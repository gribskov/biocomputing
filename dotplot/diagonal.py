"""=================================================================================================
Even though RLE provides a sparse matrix representation, it is difficult to save scores for
positions in the matching windows.  As an alternative calculate and plot one diagonal at a time

For usage examples see testing section at end

TODO: add documentation of scoring table, window, threshold, to plot and output
TODO: better legend for dots, e.g., colored sized dots

Michael Gribskov     25 June 2020
================================================================================================="""
import sys
from datetime import date
from random import choice

from bokeh.plotting import figure, show  # , output_file
from bokeh.layouts import layout
from bokeh.transform import transform
from bokeh.models import ColorBar, LinearColorMapper, ColumnDataSource, LinearAxis, Range1d, \
    BoxAnnotation
from bokeh.core.validation import silence
# noinspection PyUnresolvedReferences
from bokeh.core.validation.warnings import MISSING_RENDERERS

from sequence.score import Score
from sequence.fasta import Fasta


class Diagonal(Score, Fasta):
    """=============================================================================================

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Diagonal class constructor.
        Subclass of Score
        Delegates to Fasta via self.s1 and selfs2
        Deletates to pyplot via self.fig

        -----------------------------------------------------------------------------------------"""
        Score.__init__(self)

        self.diagonal = []
        self.threshold = 0
        self.window = 0
        self.window_max = 1
        self.window_min = 0

        # sizes of plots defined in setupBokeh()
        self.title = ''
        self.mainplot = None
        self.runplot = None
        self.scoreplot = None
        self.legend = None

        self.grid = None
        self.soff = 0  # offset so that all scores are >= zero
        self.score = None
        self.nscore = 0
        self.run = None
        self.nrun = 0

        self.rscore = None
        self.nrscore = 0
        self.rrun = None
        self.nrrun = 0

        # sequences, s1 is horizontal, s2 is vertical
        self.s1 = Fasta()
        self.s2 = Fasta()
        self.i1 = None  # integer array representation of sequences
        self.i2 = None
        self.l1 = 0
        self.l2 = 0

    def setupCalculation(self, seq1, seq2, window=5, threshold=3, resetstat=True):
        """-----------------------------------------------------------------------------------------
        Load the sequences and do some basic setup for score calculations. Sequences
        are passed as Fasta object to make it easier to use multi fasta files.

        :param seq1: Fasta object
        :param seq2: Fasta object
        :param window: int, length of window for calculation
        :param threshold: float, minimum score in window to plot
        :param resetstat: boolean, if False, reset score and run counts to zero
        :return: True
        -----------------------------------------------------------------------------------------"""
        # sequence setup
        self.s1 = seq1
        self.s2 = seq2
        self.l1, self.l2 = self.seqToInt()
        if self.l1 < self.l2:
            # shorter sequence is always s2
            self.s1, self.s2 = self.s2, self.s1
            self.l1, self.l2 = self.l2, self.l1
            self.i1, self.i2 = self.i2, self.i1

        self.diagonal = [0 for _ in range(min(self.l1, self.l2))]
        self.window = window
        self.threshold = threshold

        # stat() histograms.  nrun is always positive and in the range 0 .. window
        # maximum and minimum scores would be self.max*window amd self.max.window
        # max and min are inherited from Score
        # self.run = [0 for _ in range(window + 1)]
        if resetstat:
            score_range = int(self.max * window) - int(self.min * window) + 1
            self.soff = int(self.min * window)
            self.run = [0 for _ in range(min(self.l1, self.l2) - window + 2)]
            self.rrun = [0 for _ in range(min(self.l1, self.l2) - window + 2)]
            self.score = [0 for _ in range(score_range)]
            self.rscore = [0 for _ in range(score_range)]
            self.nscore = 0
            self.nrscore = 0
            self.nrun = 0
            self.nrrun = 0

        # score scaling
        if seq1.isACGT() and seq2.isACGT():
            # DNA sequences, range is max and min scores from table times window size
            self.window_max = window * self.max
            self.window_min = max(self.threshold - 1, 0)
        else:
            factor = 3.0
            self.window_max = factor * self.max
            self.window_min = max(self.threshold - 1, self.min)

        return True

    def setupBokeh(self):
        """-----------------------------------------------------------------------------------------
        SEt up four plot in 2 x 2 grid, but with differing sizes
            mainplot is the dotplot itself, upper right
            legend shows the colorbar legend
            scoreplot shows the window score distribution
            runploot shows the log of the run length distribution

        :return: True
        -----------------------------------------------------------------------------------------"""

        # turn off MISSING_RENDERERS warning caused by plotting colorbars in empty plot
        silence(MISSING_RENDERERS, True)

        if self.title:
            titlestr = self.title
        else:
            now = date.today()
            titlestr = 'Dotplot of {} and {} - {}'.format(self.s1.id, self.s2.id, now)

        xlabel = '\n'.join([self.s1.id, self.s1.doc])
        ylabel = '\n'.join([self.s2.doc, self.s2.id])

        # account for sequence length difference
        xlen = 800
        ylen = xlen * self.l2 / self.l1

        self.mainplot = figure(title=titlestr, x_axis_label=xlabel, y_axis_label=ylabel,
                               height=int(ylen), width=int(xlen), align='center')
        self.legend = figure(height=int(ylen), width=200)
        self.scoreplot = figure(height=300, width=500)
        self.runplot = figure(height=300, width=500, y_axis_type='log')
        self.grid = layout([[self.mainplot, self.legend], [self.scoreplot, self.runplot]])

        return True

    def setupPalette(self, base, levels, color_reverse):
        """-----------------------------------------------------------------------------------------
        Colormaps are used in multiple methods so this utility provides a unified safe method for
        setup.  Bokeh handles colormaps a little differently than other plotting programs

        :param base: string, e.g. Greys, Blues, Reds, Viridis, etc
        :param levels: int, usually 0-9 or 256
        :param color_reverse: boolean, if True highest color is dark
        :return:
        -----------------------------------------------------------------------------------------"""
        from bokeh.palettes import all_palettes

        # the defaults are here instead of in definition so that they never change
        default_base = 'Greys'
        default_level = 256
        default_reverse = True

        try:
            palette = all_palettes[base][levels]
        except (KeyError, IndexError) as error:
            # if lookup fails, use default
            palette = all_palettes[default_base][default_level]
            color_reverse = default_reverse
            sys.stderr.write('Diagonal::setupPalettes - {}, color {} levels {} is undefined.\n'.
                             format(error, base, levels))
            sys.stderr.write('\tUsing default {}{}\n'.format(default_base, default_level))

        if color_reverse:
            palette = palette[::-1]

        return palette

    def seqToInt(self):
        """-----------------------------------------------------------------------------------------
        Convert sequence strings to an integer arrays and stores in object.  An integer array is
        more convenient for direct lookups in the scoring table than a string

        :return: int, int length of sequence lists
        -----------------------------------------------------------------------------------------"""
        a2i = self.a2i

        self.i1 = [a2i[c] for c in self.s1.seq]
        self.i2 = [a2i[c] for c in self.s2.seq]

        return len(self.i1), len(self.i2)

    def rle2coord(self):
        """-----------------------------------------------------------------------------------------
        Return a list of beginning and ending positions of each run.  List is a list of four
        coordinates for each run [s1begin, s1end, s2begin, s2end]

        :return:
        -----------------------------------------------------------------------------------------"""
        coord = []
        l2 = self.l2

        for diag in range(len(self.diagonal)):

            for offset, length in self.diagonal[diag]:
                end1 = max(diag - l2 + 1, 0) + offset
                end2 = max(l2 - diag - 1, 0) + offset
                beg1 = end1 - length + 1
                beg2 = end2 - length + 1
                coord.append([beg1, end1, beg2, end2])

        return coord

    def diagLenBegin(self, diag):
        """-----------------------------------------------------------------------------------------
        Calculates the length of diagonal diag and the beginning position of the diagonal in
        each sequence

        :param diag: int, diagonal number
        :return: int (diagonal length), int (seq1 begin), int (seq2 begin)
        -----------------------------------------------------------------------------------------"""
        pos1 = max(diag - self.l2 + 1, 0)
        pos2 = max(self.l2 - diag - 1, 0)
        diaglen = min(self.l1 - pos1, self.l2 - pos2)

        return diaglen, pos1, pos2

    def diagonalScore(self, d):
        """-----------------------------------------------------------------------------------------
        Calculate the moving window sum of comparison score along one diagonal and store in the
        object.

        :param d: int, diagonal number
        :return: list, scores along diagonal
        -----------------------------------------------------------------------------------------"""
        diaglen, pos1, pos2 = self.diagLenBegin(d)

        i1 = self.i1
        i2 = self.i2

        window = self.window
        cmp = self.table
        diagonal = self.diagonal

        nmatch = 0
        old1 = pos1
        old2 = pos2

        if diaglen < window:
            # skip   diagonals shorter than window length
            return nmatch

        diagonal[:] = map(lambda i: 0, diagonal)  # lambda much faster to set all values
        # to zero
        score = 0

        # first window
        for offset in range(window):
            score += cmp[i1[pos1]][i2[pos2]]
            pos1 += 1
            pos2 += 1

        dpos = 0
        diagonal[dpos] = score

        # rest of diagonal
        for offset in range(window, diaglen):
            score -= cmp[i1[old1]][i2[old2]]
            score += cmp[i1[pos1]][i2[pos2]]

            dpos += 1
            diagonal[dpos] = score

            old1 += 1
            old2 += 1
            pos1 += 1
            pos2 += 1

        return diagonal

    def random(self, n=10000):
        """-----------------------------------------------------------------------------------------
        Calculate random score distribution using current scoring table, window, and threshold.
        Use stat() to get distributions and run lengths

        :param n: int, number of windows to calculate
        :return: list of n scores
        -----------------------------------------------------------------------------------------"""
        window = self.window
        cmp = self.table
        i1 = self.i1
        i2 = self.i2

        dist = [0 for _ in range(n - window)]
        win = [0 for _ in range(window)]

        wsum = 0
        for i in range(window):
            a = choice(i1)
            b = choice(i2)
            score = cmp[a][b]
            win[i] = score
            wsum += score

        newpos = 0
        pos = 0
        for i in range(n - window):
            dist[pos] = wsum
            wsum -= win[newpos]
            a = choice(i1)
            b = choice(i2)
            score = cmp[a][b]

            wsum += score
            win[newpos] = score

            newpos = (newpos + 1) % window
            pos += 1

        return dist

    def stat(self, diaglen, random=False, rdiagonal=None):
        """-----------------------------------------------------------------------------------------
        Save histograms of  scores and run lengths for the current diagonal. For random simulation
        create one long diagonal using the random() method.

        :param diaglen: int, length of diagonal
        :param rdiagonal: int, length of diagonal
        :param random: boolean, if True, stats are stored in random
        :return: int, int, number of runs, number of scores
        -----------------------------------------------------------------------------------------"""
        if diaglen < self.window:
            return self.nrun, self.nscore

        window = self.window
        threshold = self.threshold
        soff = self.soff

        if random:
            score = self.rscore
            nscore = self.nrscore
            run = self.rrun
            nrun = self.nrrun
            diagonal = rdiagonal
        else:
            score = self.score
            nscore = self.nscore
            run = self.run
            nrun = self.nrun
            diagonal = self.diagonal

        runlen = 0
        for offset in range(diaglen - window):
            score[int(diagonal[offset]) - soff] += 1
            nscore += 1
            if diagonal[offset] >= threshold:
                runlen += 1

            else:
                run[runlen] += 1
                runlen = 0
                nrun += 1

        if runlen:
            # print(runlen)
            run[runlen] += 1
            nrun += 1

        if random:
            self.nrscore = nscore
            self.nrrun = nrun
        else:
            self.nscore = nscore
            self.nrun = nrun

        return nrun, nscore

    def statPlot(self, write_cumulative=True):
        """-----------------------------------------------------------------------------------------
        Plot the score distribution and the log of the run lengths (with +1 prior)
        TODO: hover to get exact values?

        :param write_cumulative: boolean, if true (default) cumulative distribution is written to
        STDOUT
        :return: True
        -----------------------------------------------------------------------------------------"""
        scoreplot = self.scoreplot
        runplot = self.runplot

        # ------------------------
        # random simulation - approximately same size as real comparison
        # ------------------------
        nr = int((self.l1 - self.window) * (self.l2 - self.window))
        random = self.random(n=nr)
        self.stat(nr, random=True, rdiagonal=random)
        self.density(self.rscore, self.nrscore)

        # ------------------------
        # score distribution
        # ------------------------
        score = self.score
        soff = self.soff
        nscore = self.nscore

        # find the first and last non-zero index in the score frequency data
        # make this the min and max for the x axis, adjust x values for the score offset
        scoremin, scoremax = self.scoreMinMax(score)
        maxp = self.density(score, nscore)
        cumulative = self.cumulative(score, 1)
        x = [i + soff for i in range(scoremin, scoremax + 1)]

        if write_cumulative:
            sys.stdout.write('Cumulative score distribution\n\tScore\t\tProb\t\tDensity\n')

            for i in range(len(x)):
                j = i + scoremin
                sys.stdout.write('\t{}\t\t\t{:.4f}\t\t{:.4f}\n'.format(x[i], cumulative[j],
                                                                       score[j]))

        # observed score density (blue)
        scoreplot.vbar(x=x, top=score[scoremin:scoremax + 1], width=0.8, color='#0000ff',
                       line_color='black', alpha=0.5, bottom=0.0)
        scoreplot.y_range = Range1d(0.0, maxp * 1.1)

        # cumulative distribution, blue line on right side axis
        scoreplot.extra_y_ranges = {"cumulative": Range1d(start=0.0, end=1.0)}
        axis2 = LinearAxis(y_range_name="cumulative")
        axis2.ticker.num_minor_ticks = 10
        scoreplot.add_layout(axis2, 'right')
        scoreplot.line(x=x, y=cumulative[scoremin:scoremax + 1], y_range_name='cumulative',
                       line_width=2, color='#1122cc')

        # shaded box showing 95% level
        box = BoxAnnotation(bottom=0.95, top=1.0, y_range_name='cumulative',
                            fill_color='#FFBBBB', line_width=3, line_dash='dashed')
        scoreplot.add_layout(box)

        # random score density (blue)
        scoreplot.vbar(x=x, top=self.rscore[scoremin:scoremax + 1], width=0.8, color='#ff0000',
                       line_color='black', alpha=0.5, bottom=0.0)

        # alternative plot of random density as a line, i didn't like it
        # scoreplot.line(x=x, y=self.rscore[scoremin:scoremax + 1], width=0.8, color='#FF3333',
        #                line_width=2)

        # ------------------------
        # run distribution
        # ------------------------
        run = self.run
        rrun = self.rrun
        maxrun = 0
        for i in range(max(len(run), len(rrun))):
            if run[i] or rrun[i] > 0:
                maxrun = i
            # + 1 to prevent log errors
            run[i] += 1
            rrun[i] += 1

        minrun = 1
        x = [i for i in range(minrun, maxrun + 1)]
        # observed and simulated run lengths,  need bottom=1 because of log axis

        runplot.vbar(x=x, top=run[minrun:maxrun + 1], width=0.8,
                     color='#0000ff', line_color='black', alpha=0.5,
                     line_width=0.5, bottom=1.0)

        runplot.vbar(x=x, top=rrun[minrun:maxrun + 1], width=0.8,
                     color='#ff0000', line_color='black', alpha=0.5,
                     line_width=0.5, bottom=1.0)

        return True

    def drawDot(self, rev=False, cbase='Greys', clevel=256, crev=True, width=True,
                color=True, alpha=0.5):
        """-----------------------------------------------------------------------------------------
        Draw all diagonals as dots.  the score is reported in the first position of the window so
        the position of the dot must be offset to lie in the middle of the window.  Zero origin
        coordinates must also be incremented by 1

        :return: True
        -----------------------------------------------------------------------------------------"""

        window = self.window
        halfwindow = (window - 1) / 2.0
        threshold = self.threshold
        l2 = self.l2

        cmap = LinearColorMapper(palette=self.setupPalette(cbase, clevel, crev),
                                 low=self.window_min, high=self.window_max)
        col_max = window

        # reversed plot: invert the direction and change color
        # markers are slightly smaller on reversed plot so they show ujp when superimposed
        yinc = 1
        if rev:
            yinc = -1
            # size_max /= 1.25

        size_min = 2
        size_max = 10

        score_min = self.window_min
        score_max = self.window_max

        try:
            wscale = (size_max - size_min) / (score_max - score_min)
            woff = size_min - score_min * wscale
        except ZeroDivisionError:
            wscale = 1
            woff = size_min

        for d in range(self.l1 + self.l2 - 1):
            dscore = self.diagonalScore(d)
            diaglen, xpos, ypos = self.diagLenBegin(d)
            self.stat(diaglen)

            xpos += halfwindow
            if rev:
                ypos = l2 - ypos - halfwindow - 1
            else:
                ypos += halfwindow

            # plot one diagonal at a time, less memory?
            dot = {'x': [], 'y': [], 'size': [], 'score': []}
            source = ColumnDataSource(data=dot)

            for pos in range(diaglen - window + 1):
                if dscore[pos] >= threshold:
                    size = size_min
                    if width:
                        size = min(dscore[pos] * wscale + woff, size_max)

                    col = col_max
                    if color:
                        col = dscore[pos]

                    dot['x'].append(xpos)
                    dot['y'].append(ypos)
                    dot['size'].append(size)
                    dot['score'].append(col)

                xpos += 1
                ypos += yinc

            self.mainplot.circle(x='x', y='y', source=source,
                                 size='size',
                                 line_color=transform('score', cmap), line_alpha=alpha,
                                 fill_color=transform('score', cmap), fill_alpha=alpha)

        # color bar is in a separate window, self.legend, so it doesn't disturb the
        # aspect ratio
        color_bar = ColorBar(color_mapper=cmap, label_standoff=3, bar_line_color='black',
                             scale_alpha=alpha, width=20, margin=0, location=(0, 0),
                             major_tick_in=20, major_tick_out=5, major_tick_line_color='black')

        self.legend.add_layout(color_bar, 'left')

        return True

    def drawSegment(self, rev=False, cbase='Greys', clevel=256, crev=True, width=True,
                    color=True, dither=0.5, alpha=0.5):
        """-----------------------------------------------------------------------------------------
        Draw all diagonals as dots.  the score is reported in the first position of the window so
        the position of the dot must be offset to lie in the middle of the window.  Zero origin
        coordinates must also be incremented by 1

        :return: True
        -----------------------------------------------------------------------------------------"""

        window = self.window
        halfwindow = (window - 1) / 2.0
        threshold = self.threshold
        l2 = self.l2

        cmap = LinearColorMapper(palette=self.setupPalette(cbase, clevel, crev),
                                 low=self.window_min, high=self.window_max)
        col_max = window

        # reversed plot: invert the direction
        # dither sets up the segments to run from halfawy between character positions
        yinc = 1
        ydither = [-dither, + dither]
        if rev:
            yinc = -1
            ydither = [dither, -dither]

        size_min = 1
        size_max = 8

        score_min = self.window_min
        score_max = self.window_max

        try:
            wscale = (size_max - size_min) / (score_max - score_min)
            woff = size_min - score_min * wscale
        except ZeroDivisionError:
            wscale = 1
            woff = size_min

        for d in range(self.l1 + self.l2 - 1):
            dscore = self.diagonalScore(d)
            diaglen, xpos, ypos = self.diagLenBegin(d)
            self.stat(diaglen)

            xpos += halfwindow
            if rev:
                ypos = l2 - ypos - halfwindow - 1
            else:
                ypos += halfwindow

            # plot one diagonal at a time, less memory?
            segment = {'x0': [], 'y0': [], 'x1': [], 'y1': [], 'size': [], 'score': []}
            source = ColumnDataSource(data=segment)

            for pos in range(diaglen - window + 1):
                if dscore[pos] >= threshold:
                    size = size_min
                    if width:
                        size = min(dscore[pos] * wscale + woff, size_max)

                    col = col_max
                    if color:
                        col = dscore[pos]

                    segment['x0'].append(xpos - dither)
                    segment['y0'].append(ypos + ydither[0])
                    segment['x1'].append(xpos + dither)
                    segment['y1'].append(ypos + ydither[1])
                    segment['size'].append(size)
                    segment['score'].append(col)

                xpos += 1
                ypos += yinc

            self.mainplot.segment(x0='x0', x1='x1', y0='y0', y1='y1', alpha=alpha, source=source,
                                  line_width='size', line_color=transform('score', cmap))

        # color bar is in a separate window, self.legend, so it doesn't disturb the
        # aspect ratio
        color_bar = ColorBar(color_mapper=cmap, label_standoff=3, bar_line_color='black',
                             scale_alpha=alpha, width=20, margin=0, location=(0, 0),
                             major_tick_in=20, major_tick_out=5, major_tick_line_color='black')

        self.legend.add_layout(color_bar, 'right')

        return True

    def show(self, *args, **kwargs):
        """-----------------------------------------------------------------------------------------
        Delegate to plt.show().  Makes syntax a little easier in application since the object is
        used instead of the matplotlib class (application doesn't have to know matplotlib).

        :param args:
        :param kwargs:
        :return: True
        -----------------------------------------------------------------------------------------"""
        show(self.grid, *args, **kwargs)

        return True

    @staticmethod
    def cumulative(score, total):
        """-----------------------------------------------------------------------------------------
        Return cumulative score probability distribution as a list

        :param score: list
        :param total: int, number of observations
        :return: list
        -----------------------------------------------------------------------------------------"""
        cumulative = []
        wsum = 0
        for i in range(len(score)):
            wsum += score[i] / total
            cumulative.append(wsum)

        return cumulative

    @staticmethod
    def density(score, total):
        """-----------------------------------------------------------------------------------------

        :param score:
        :param total:
        :return:
        -----------------------------------------------------------------------------------------"""
        maxp = 0.0
        for i in range(len(score)):
            score[i] /= total
            maxp = max(maxp, score[i])

        return maxp

    @staticmethod
    def scoreMinMax(score):
        """-----------------------------------------------------------------------------------------
        Returns the first and last non-zero positions in a list of scores.  Use to get ranges for
        score histograms

        :param score: list
        :return: int, int
        -----------------------------------------------------------------------------------------"""
        scoremin = None
        scoremax = None
        for i in range(len(score)):
            if score[i] > 0:
                if scoremin is None:
                    scoremin = i
                scoremax = i

        return scoremin, scoremax


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    match = Diagonal()

    fasta1 = Fasta(filename=sys.argv[1])
    fasta1.read()

    fasta2 = Fasta()
    fasta2.seq = fasta1.seq
    fasta2.id = 'seq 2'
    fasta2.doc = 'Sequence 2'

    # select list of tests to run, one plot in browser for each
    # tests = [0,1,2,3]
    # tests = range(4,10)
    tests = [8, 9]

    for test in tests:

        if test == 0:
            w = 1
            t = 1
            match.title = 'test {} - forward dotplot size cueing only, w={} t={}'. \
                format(test, w, t)
            match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
            match.setupBokeh()
            match.drawDot(width=True, color=False, alpha=1)
            match.statPlot()

        elif test == 1:
            w = 3
            t = 1
            match.title = 'test {} - forward dotplot size cueing only, w={} t={}'. \
                format(test, w, t)
            match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
            match.setupBokeh()
            match.drawDot(width=True, color=False, alpha=1)
            match.statPlot()

        elif test == 2:
            w = 12
            t = 8
            match.title = 'test {} - different lengths, forward dots, no cueing, w={} t={}'. \
                format(test, 1, 1)
            fasta1.seq = fasta1.seq[:200]
            fasta2.seq = fasta1.seq[:50]
            fasta2.doc = 'bases 1 to 50'
            fasta2.doc = 'bases 1 to 200'
            match.setupCalculation(fasta1, fasta2, window=1, threshold=1)
            match.setupBokeh()
            match.drawDot(width=False, color=False)
            match.statPlot()

        elif test == 3:
            w = 12
            t = 8
            match.title = 'test {} - different lengths, short on  axis, w={} t={}'. \
                format(test, 1, 1)
            fasta2.seq = fasta1.seq[:200]
            fasta1.seq = fasta1.seq[:50]
            fasta1.doc = 'bases 1 to 50'
            fasta2.doc = 'bases 1 to 200'
            match.setupCalculation(fasta2, fasta1, window=1, threshold=1)
            match.setupBokeh()
            match.drawDot(width=False, color=False)
            match.statPlot()

        elif test == 4:
            w = 10
            t = 8
            match.title = 'test {} - forward dotplot color cueing only, w={} t={}'. \
                format(test, w, t)
            match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
            match.setupBokeh()
            match.drawDot(cbase='Viridis', clevel=256, crev=True, width=False, color=True, alpha=1)
            match.statPlot()

        elif test == 5:
            w = 12
            t = 8
            match.title = 'test {} - forward/reverse dotplot, w={} t={}'. \
                format(test, w, t)
            match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
            match.setupBokeh()
            match.drawDot(cbase='Blues', clevel=256, crev=True, width=True)
            fasta2.seq = fasta1.reverseComplement()
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t, resetstat=False)
            match.drawDot(color=True, cbase='Reds', clevel=256, crev=True, rev=True, width=True)
            match.statPlot()

        elif test == 6:
            w = 12
            t = 8
            match.title = 'test {} - forward/reverse lines, width only, w={} t={}'. \
                format(test, w, t)
            fasta2.seq = fasta1.reverseComplement()
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh()
            match.drawSegment(cbase='Blues', clevel=256, crev=True, rev=True, width=False,
                              color=False)
            match.statPlot()

        elif test == 7:
            w = 12
            t = 8
            match.title = 'test {} - forward/reverse lines, width and color, w={} t={}'. \
                format(test, w, t)
            match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
            match.setupBokeh()
            match.drawSegment(color=True, cbase='Blues', clevel=256, crev=True, width=True)
            fasta2.seq = fasta1.reverseComplement()
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t, resetstat=False)
            match.drawSegment(color=True, cbase='Reds', clevel=256, crev=True, rev=True, width=True,
                              alpha=0.8)
            match.statPlot()

        elif test == 8:
            w = 12
            t = 9
            match.title = 'test {} - protein dot with blosum, w={} t={}'. \
                format(test, w, t)
            match.readNCBI('table/BLOSUM62.matrix')

            fasta = Fasta(filename='../data/globin.pfa')
            fasta.read()
            fasta1 = fasta.copy()
            fasta.read()
            fasta2 = fasta.copy()

            match.setupCalculation(fasta2, fasta1, window=w, threshold=t)
            match.setupBokeh()
            match.drawDot(cbase='Viridis', clevel=256, crev=True)
            match.statPlot()

        elif test == 9:
            w = 18
            t = 12
            match.title = 'test {} - protein line with blosum, w={} t={}'. \
                format(test, w, t)
            match.readNCBI('table/BLOSUM62.matrix')

            fasta = Fasta(filename='../data/globin.pfa')
            fasta.read()
            fasta1 = fasta.copy()
            fasta.read()
            fasta2 = fasta.copy()

            match.setupCalculation(fasta2, fasta1, window=w, threshold=t)
            match.setupBokeh()
            match.drawSegment(cbase='Viridis', clevel=256, crev=True)
            match.statPlot()
        match.show()

exit(0)

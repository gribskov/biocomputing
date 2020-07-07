"""=================================================================================================
Even though RLE provides a sparse matrix representation, it is difficult to save scores for
positions in the matching windows.  As an alternative calculate and plot one diagonal at a time

TODO: need to deal with non-identity scoring matrices

Michael Gribskov     25 June 2020
================================================================================================="""
import sys
from datetime import date

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import grid, gridplot, layout
from bokeh.transform import transform, linear_cmap
from bokeh.models import ColorBar, LinearColorMapper, ColumnDataSource
from bokeh.core.validation import silence
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
        self.ax = None
        self.title = ''
        self.threshold = 0
        self.window = 0

        # sizes of histograms defined in setup
        self.mainplot = None
        self.runplot = None
        self.scoreplot = None
        self.grid = None
        self.soff = 0  # offset so thatall scores are >= zero
        self.nrun = 0
        self.nscore = 0

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
        :return: True
        -----------------------------------------------------------------------------------------"""
        # sequence setup
        self.s1 = seq1
        self.s2 = seq2
        self.l1, self.l2 = self.seqToInt()
        self.diagonal = [0 for _ in range(min(self.l1, self.l2))]
        self.window = window
        self.threshold = threshold

        # stat() histograms.  nrun is always positive and in the range 0 .. window
        # maximum and minimum scores would be self.max*window amd self.max.window
        # max and min are inherited from Score
        # self.run = [0 for _ in range(window + 1)]
        if resetstat:
            self.run = [0 for _ in range(min(self.l1, self.l2) - window + 2)]
            score_range = int(self.max * window) - int(self.min * window) + 1
            self.soff = int(self.min * window)
            self.score = [0 for _ in range(score_range)]

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

        # turn of MISSING_RENDERERS warning caused by plotting colorbars in empty plot
        silence(MISSING_RENDERERS, True)

        if self.title:
            titlestr = self.title
        else:
            now = date.today()
            titlestr = 'Dotplot of {} and {} - {}'.format(self.s1.id, self.s2.id, now)

        xlabel = '\n'.join([self.s1.id, self.s1.doc])
        ylabel = '\n'.join([self.s2.doc, self.s2.id])

        # TODO need to set up the plot size to account for sequnce lengths

        self.mainplot = figure(title=titlestr, x_axis_label=xlabel, y_axis_label=ylabel,
                               height=800, width=800)
        self.legend = figure(height=800, width=200)

        # background_fill_color='#bbbbff',
        # aspect_ratio='auto', x_range=(1, 10))

        self.scoreplot = figure(height=300, width=500)
        self.runplot = figure(height=300, width=500, y_axis_type='log')
        # self.runplot = figure(height=300, width=500)

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

        default_base = 'Greys'
        default_level = 256
        default_reverse = True

        try:
            palette = all_palettes[base][levels]
        except (KeyError, IndexError) as error:
            palette = all_palettes[default_base][default_level]
            sys.stderr.write(
                'Diagonal::setupPalettes - color {} levels {} is undefined.\n'. \
                    format(base, levels))
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

        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = self.l1
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

    def stat(self, diaglen):
        """-----------------------------------------------------------------------------------------
        Save histograms of  scores and run lengths for the current diagonal

        Poaram diagonal: int, length of diagonal
        :return: int, int, number of runs, number of scores
        -----------------------------------------------------------------------------------------"""
        if diaglen < self.window:
            return self.nrun, self.nscore

        window = self.window
        threshold = self.threshold
        score = self.score
        nscore = self.nscore
        soff = self.soff
        run = self.run
        nrun = self.nrun
        diagonal = self.diagonal

        runlen = 0
        for offset in range(diaglen - window + 1):
            score[int(diagonal[offset]) + soff] += 1
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

        return nrun, nscore

    def statPlot(self):
        """-----------------------------------------------------------------------------------------
        Plot the score distribution and the log of the run lengths (with +1 prior)

        :return:
        -----------------------------------------------------------------------------------------"""
        scoreplot = self.scoreplot
        runplot = self.runplot

        # score distribution
        x = [i for i in range(self.window + 1)]
        scoreplot.vbar(x=x, top=self.score, width=0.8, color='#AAAAFF', line_color='black')


        # run distribution
        run = self.run
        maxrun = 0
        for i in range(len(run)):
            if run[i] > 0:
                maxrun = i
            # run[i] = log10(run[i]+1)
            run[i] += 1

        x = [i for i in range(1, maxrun + 1)]
        # need bottom=1 becasue of log axis
        runplot.vbar(x=x, top=run[1:maxrun + 1], width=0.8, color='#FF5555',
                     line_color='black', line_width=0.5, bottom=1)

        return True

    def drawDot(self, rev=False, cbase='Greys', clevel=256, crev='True', width=False,
                color=False, alpha=0.5):
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
                                 low=max(threshold-1,0), high=window)

        # reversed plot: invert the direction and change color
        # markers are slightly smaller on reversed plot so they show ujp when superimposed
        yinc = 1
        if rev:
            yinc = -1
            # size_max /= 1.25

        size_min = 2
        size_max = 8

        try:
            wscale = (size_max-size_min) / (window - threshold)
            woff = size_min - threshold * wscale
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
                    f = dscore[pos] * wscale + woff
                    size = size_min
                    if width:
                        size = f
                    if not color:
                        f = 1.0
                    dot['x'].append(xpos)
                    dot['y'].append(ypos)
                    dot['size'].append(size)
                    dot['score'].append(dscore[pos])

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

    def drawSegment(self, rev=False, cbase='Greys', clevel=256, crev='True', width=False,
                    color=False, dither=0.5, alpha=0.5):
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
                                 low=max(threshold-1,0), high=window)

        # reversed plot: invert the direction
        # dither sets up the segments to run from halfawy between character positions
        yinc = 1
        ydither = [-dither, + dither]
        if rev:
            yinc = -1
            ydither = [dither, -dither]

        size_min = 2
        size_max = 6

        try:
            wscale = (size_max - size_min) / (window - threshold)
            woff = size_min - threshold * wscale
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
                    f = dscore[pos] * wscale + woff
                    size = size_min
                    if width:
                        size = f
                    if not color:
                        f = dscore[pos] = window

                    segment['x0'].append(xpos - dither)
                    segment['y0'].append(ypos + ydither[0])
                    segment['x1'].append(xpos + dither)
                    segment['y1'].append(ypos + ydither[1])
                    segment['size'].append(size)
                    segment['score'].append(dscore[pos])
                    inwindow = True

                xpos += 1
                ypos += yinc

            self.mainplot.segment(x0='x0', x1='x1', y0='y0', y1='y1', alpha=alpha, source=source,
                                  line_width='size', line_color=transform('score', cmap))

        # color bar is in a separate window, self.legend, so it doesn't disturb the
        # aspect ratio
        color_bar = ColorBar(color_mapper=cmap, label_standoff=3, bar_line_color='black',
                             scale_alpha = alpha, width=20, margin=0, location=(0, 0),
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
        show(self.grid)


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    match = Diagonal()

    fasta1 = Fasta(filename=sys.argv[1])
    fasta1.read()

    fasta2 = Fasta()
    fasta2.id = 'seq2'
    fasta2.doc = ' bases 1:50'
    fasta2.seq = fasta1.seq[:50]

    fasta1.seq = fasta1.seq[:200]

    w = 12
    t = 7

    test = 4

    if test == 0:
        match.title = 'test 0 - forward dotplot'
        match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
        match.setupBokeh()
        match.drawDot(width=True)
        match.statPlot()

    elif test == 1:
        match.title = 'test1 - forward/reverse dotplot'
        match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
        match.setupBokeh()
        match.drawDot(cbase='Blues', clevel=256, crev=True, width=True)
        fasta2.seq = fasta1.reverseComplement()
        match.setupCalculation(fasta1, fasta2, window=w, threshold=t, resetstat=False)
        match.drawDot(color=True, cbase='Reds', clevel=256, crev=True, rev=True, width=True)
        match.statPlot()

    elif test == 2:
        match.title = 'test 2 - reverse dotplot with lines'
        fasta2.seq = fasta1.reverseComplement()
        match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
        match.setupBokeh()
        match.drawSegment(cbase='Blues', clevel=256, crev=True, rev=True, width=False,
                          color=False)
        match.statPlot()

    elif test == 3:
        match.title = 'test 3 - forward/reverse with lines'
        match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
        match.setupBokeh()
        match.drawSegment(color=True, cbase='Blues', clevel=256, crev=True, width=True)
        fasta2.seq = fasta1.reverseComplement()
        match.setupCalculation(fasta1, fasta2, window=w, threshold=t, resetstat=False)
        match.drawSegment(color=True, cbase='Reds', clevel=256, crev=True, rev=True, width=True,
                          alpha=0.8)
        match.statPlot()

    elif test == 4:
        match.title = 'test4 - self dotplot without window'
        match.setupCalculation(fasta1, fasta1, window=1, threshold=1)
        match.setupBokeh()
        match.drawDot(width=False, color=False)
        match.statPlot()


    match.show()

    exit(0)

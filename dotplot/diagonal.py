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
    BoxAnnotation, Scatter
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
        Delegates to Fasta via self.s1 and self.s2
        Delegates to pyplot via self.fig

        diagonal: one diagonal of scores
        yinc: direction of y axis, 1 or -1 means forward or reverse respectively
        window: length of window for calculation
        threshold: minimum value for window to be plotted

        -----------------------------------------------------------------------------------------"""
        Score.__init__(self)

        self.diagonal = []
        self.single = False
        self.yinc = 1
        self.threshold = 0
        self.window = 0
        self.nscore = 0
        self.nrun = 0

        self.frame = {}  # data frames
        self.function = {}  # functions for populating data frames

        # Plotting variables
        # sizes of panels are defined in setupBokeh()
        self.title = ''
        self.figure = {}
        self.grid = None
        self.palette = None
        self.cmap = None
        self.alpha = 0.5
        self.mindotsize = 2
        self.maxdotsize = 10

        # sequences, s1 is horizontal, s2 is vertical
        self.s1 = Fasta()
        self.s2 = Fasta()
        self.i1 = None  # integer array representation of sequences
        self.i2 = None
        self.l1 = 0
        self.l2 = 0
        self.seqreverse = False  # only applies to s2

    def setupCalculation(self, seq1, seq2, window=5, threshold=3, resetstat=True):
        """-----------------------------------------------------------------------------------------
        Load the sequences and do some basic setup for score calculations. Sequences
        are passed as Fasta object to make it easier to use multi fasta files.

        :param seq1: hash of id, doc, seq, format
        :param seq2:
        :param window: int, length of window for calculation
        :param threshold: float, minimum score in window to plot
        :param resetstat: boolean, if False, reset score and run counts to zero
        :return: True
        -----------------------------------------------------------------------------------------"""
        # sequence setup
        self.s1 = seq1
        self.s2 = seq2
        self.l1 = len(seq1['seq'])
        self.l2 = len(seq2['seq'])

        # move shorter sequence to s2 if necessary
        if self.l1 < self.l2:
            # shorter sequence is always s2
            self.s1, self.s2 = self.s2, self.s1
            self.l1, self.l2 = self.l2, self.l1

        # reverse sequence 2 if necessary

        if self.seqreverse:
            if self.s2['dir'] == 'forward':
                # reverse forward sequence
                self.s2['seq'] = Fasta.reverseComplement(self.s2['seq'])
                self.s2['dir'] = 'reverse'
            self.yinc = -1
        else:
            if self.s2['dir'] == 'reverse':
                # un-reverse previously reversed
                self.s2['seq'] = Fasta.reverseComplement(self.s2['seq'])
                self.s2['dir'] = 'forward'
            self.yinc = 1

        # setup integer array version of sequence
        self.seqToInt()

        self.diagonal = [0 for _ in range(min(self.l1, self.l2))]
        self.window = window
        self.threshold = threshold

        # stat() histograms.  nrun is always positive

        if resetstat:
            self.nscore = 0
            self.nrun = 0
            for frame in self.frame:
                self.resetFrame(frame)

        return True

    def setupBokeh(self, cbase=None, clevels=None, creverse=None):
        """-----------------------------------------------------------------------------------------
        SEt up four plot in 2 x 2 grid, but with differing sizes
            mainplot is the dotplot itself, upper right
            legend shows the colorbar legend
            scoreplot shows the window score distribution
            runploot shows the log of the run length distribution

        :param cbase: string, e.g. Greys, Blues, Reds, Viridis, etc
        :param clevels: int, usually 0-9 or 256
        :param creverse: boolean, if True highest color is dark
        :return: True
        -----------------------------------------------------------------------------------------"""

        # turn off MISSING_RENDERERS warning caused by plotting colorbars in empty plot
        silence(MISSING_RENDERERS, True)

        self.palette = self.setupPalette(cbase=cbase, clevels=clevels, creverse=creverse)

        if self.title:
            titlestr = self.title
        else:
            now = date.today()
            titlestr = 'Dotplot of {} and {} - {}'.format(self.s1['id'], self.s2['id'], now)

        xlabel = '\n'.join([self.s1['id'], self.s1['doc']])
        ylabel = '\n'.join([self.s2['doc'], self.s2['id']])

        # account for sequence length difference, ylen scaling affects main and legend panels
        xlen = 800
        ylen = xlen * self.l2 / self.l1

        # define each panel as a figure
        label = '({}, {}, score)'.format(self.s1['id'], self.s2['id'])
        TIPS = [(label, '($x{0}, $y{0}, @score)')]
        self.figure['main'] = figure(title=titlestr, x_axis_label=xlabel, y_axis_label=ylabel,
                                     height=int(ylen), width=int(xlen), align='center',
                                     tooltips=TIPS)

        self.figure['legend'] = figure(height=int(ylen), width=200)

        TIPS = [('score, number', '$x{0}, $y{0.00}')]
        self.figure['scoredist'] = figure(height=300, width=500, tooltips=TIPS)

        TIPS = [('length,count', '$x{0}, $y{0}')]
        self.figure['rundist'] = figure(height=300, width=500, y_axis_type='log', tooltips=TIPS)

        # grid layout
        self.grid = layout([[self.figure['main'], self.figure['legend']],
                            [self.figure['scoredist'], self.figure['rundist']]])

        return True

    def setupFrame(self, defs):
        """-----------------------------------------------------------------------------------------
        Setup data frames for the defined analyses with empty ndata fields.  Each def in defs
        defines
            name - name of data frame
            function - a callback function used to construct the data from a diagonal of scores
            variables - variables that will be populated

        As used here, a dataframes are stored in the object as self.frame[name]
        self.frame[name] = {function, var1: [], var2: [], var3: [], ...}

        :param defs: list, see above
        :return: int, number of frames define
        -----------------------------------------------------------------------------------------"""
        n = 0
        for defin in defs:
            n += 1
            self.frame[defin['data']] = {}
            self.function[defin['data']] = defin['fn']
            for v in defin['var']:
                self.frame[defin['data']][v] = []

        return n

    def resetFrame(self, framename):
        """-----------------------------------------------------------------------------------------
        Reset the data in one frame to empty lists.  Needed for reverse plots

        :param framename:
        :return: True
        -----------------------------------------------------------------------------------------"""
        frame = self.frame[framename]
        for var in frame:
            frame[var] = []

        return True

    def setupPalette(self, cbase, clevels, creverse):
        """-----------------------------------------------------------------------------------------
        Colormaps are used in multiple methods so this utility provides a unified safe method for
        setup.  Bokeh handles colormaps a little differently than other plotting programs

        :param cbase: string, e.g. Greys, Blues, Reds, Viridis, etc
        :param clevels: int, usually 0-9 or 256
        :param creverse: boolean, if True highest color is dark
        :return:
        -----------------------------------------------------------------------------------------"""
        from bokeh.palettes import all_palettes

        # the defaults are here instead of in definition so that they never change
        default_base = 'Greys'
        default_levels = 256
        default_reverse = True

        try:
            palette = all_palettes[cbase][clevels]
        except (KeyError, IndexError) as error:
            # if lookup fails, use default
            palette = all_palettes[default_base][default_levels]
            creverse = default_reverse
            sys.stderr.write('Diagonal::setupPalettes - {}, color {} levels {} is undefined.\n'.
                             format(error, cbase, clevels))
            sys.stderr.write('\tUsing default {}{}\n'.format(default_base, default_levels))

        if creverse:
            # reverse the orde of colors
            palette = palette[::-1]

        return palette

    def seqToInt(self):
        """-----------------------------------------------------------------------------------------
        Convert sequence strings to an integer arrays and stores in object.  An integer array is
        more convenient for direct lookups in the scoring table than a string

        :return: int, int length of sequence lists
        -----------------------------------------------------------------------------------------"""
        a2i = self.a2i

        self.i1 = [a2i[c] for c in self.s1['seq']]
        self.i2 = [a2i[c] for c in self.s2['seq']]

        return len(self.i1), len(self.i2)

    def rle2coord(self):
        """-----------------------------------------------------------------------------------------
        Return a list of beginning and ending positions of each run.  List is a list of four
        coordinates for each run [s1begin, s1end, s2begin, s2end]

        :return: 4 x int, beg1, end1, beg2, end2
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

        # if self.seqreverse:
        # pos2 = self.l2 - pos2

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

        old1 = pos1
        old2 = pos2

        if diaglen < window:
            # skip   diagonals shorter than window length
            return []

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
            # sys.stderr.write('{}\t{}\n'.format(pos1,pos2))
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
        Use stat() to get distributions and run lengths.  Use n = number of windows calculated for
        actual sequences.

        :param n: int, number of windows to calculate
        :return: list of n scores
        -----------------------------------------------------------------------------------------"""
        window = self.window
        cmp = self.table
        i1 = self.i1
        i2 = self.i2

        if n == 0:
            n = self.l1 * self.l2

        self.diagonal = [0 for _ in range(n - window)]
        dist = self.diagonal
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

    def allDiagonals(self, select):
        """-----------------------------------------------------------------------------------------
        Iterate over all diagonals and apply specified actions to each diagonal.  Each action is
        a tuple that specifies the name of the resulting data frame, and a function to process
        the diagonal. The frames are usable as Bokeh sources for plotting.

        :param select: list, names of dataframes to calculate from each diagonal
        :return: True
        -----------------------------------------------------------------------------------------"""
        frame = self.frame
        function = self.function

        for d in range(self.l1 + self.l2 - 1):
            dscore = self.diagonalScore(d)
            if not dscore:
                continue

            for data in select:
                # apply each selected function to this diagonal of scores to populate the
                # dataframes
                fxn = function[data]
                fxn(data, d)

        return True

    def windowThreshold(self, framename, d):
        """-----------------------------------------------------------------------------------------
        Callback function for allDiagonals.  Savs windows with score >= threshold in dataframe
        framename.  Works on the internally stored diagonal of scores calculated by diagonalScore()

        :param framename: string, name of a dataframe in self.frame
        :param d: int, diagonal number
        :return: True
        -----------------------------------------------------------------------------------------"""
        frame = self.frame[framename]
        dscore = self.diagonal
        window = self.window
        halfwindow = (window - 1) / 2.0
        threshold = self.threshold
        yinc = self.yinc

        diaglen, xpos, ypos = self.diagLenBegin(d)
        if diaglen < window:
            return False

        xpos += halfwindow
        if self.yinc < 0:
            ypos = self.l2 - ypos - halfwindow - 1
        else:
            ypos += halfwindow

        for pos in range(diaglen - window + 1):
            if dscore[pos] >= threshold:
                frame['x'].append(xpos)
                frame['y'].append(ypos)
                frame['score'].append(dscore[pos])

            xpos += 1
            ypos += yinc
            self.nscore += 1

        return True

    def scaleColumn(self, framename, column_source, column_dest, value, scale):
        """-----------------------------------------------------------------------------------------
        Performs a simple linear scaling on a column

        :param framename: string, a data frame in self.frame
        :param column_source: string, the column in frame to be scaled
        :param column_dest: string, name for the scaled column (in frame)
        :param value: tuple, low and high value for the input data
        :param scale: tuple, low and high value for the scaled data
        :return:
        -----------------------------------------------------------------------------------------"""
        frame = self.frame[framename]
        values = frame[column_source]
        frame[column_dest] = []
        # width = frame[column_dest]

        rangeval = value[1] - value[0]
        rangesize = scale[1] - scale[0]
        m = rangesize / rangeval

        for v in values:
            size = scale[0] + (v - value[0]) * m
            frame[column_dest].append(size)

        return

    def histogramScore(self, scoreframe, d):
        """-----------------------------------------------------------------------------------------
        Callback function for allDiagonals. Creates data frames with the score distribution. Works
        on the internally stored diagonal of scores calculated by diagonalScore()

        :param scoreframe: string, name of dataframe in self.frame
        :param d: int, diagonal number
        :return: int, number of values in columns of dataframe
        -----------------------------------------------------------------------------------------"""
        scoreframe = self.frame[scoreframe]
        diagonal = self.diagonal
        window = self.window

        if self.single:
            diaglen = len(diagonal)
        else:
            diaglen, xpos, ypos = self.diagLenBegin(d)
            diaglen -= window - 1

        nscore = 0
        score = {}
        for s in diagonal[:diaglen]:
            try:
                score[s] += 1
            except KeyError:
                score[s] = 1

            nscore += 1

        # insert into data frame, the dateframe is randomly ordered
        for s in score:
            try:
                i = scoreframe['score'].index(s)
                scoreframe['count'][i] += score[s]

            except ValueError:
                scoreframe['score'].append(s)
                scoreframe['count'].append(score[s])

        return len(scoreframe['score'])

    def histogramRun(self, runframe, d):
        """-----------------------------------------------------------------------------------------
        Callback function for allDiagonals. Create a dataframe with the run length distribution,
        apply the threshold stored in self.threshold.
        Works on the internally stored diagonal of scores calculated by diagonalScore()

        :param runframe: string, name of dataframe in self.frame
        :param d: int, diagonal number
        :return: int, number of values in columns of dataframe
        -----------------------------------------------------------------------------------------"""
        runframe = self.frame[runframe]
        diagonal = self.diagonal
        window = self.window
        threshold = self.threshold

        if self.single:
            diaglen = len(diagonal)
        else:
            diaglen, xpos, ypos = self.diagLenBegin(d)
            diaglen -= window - 1

        run = {}
        nrun = 0
        runlen = 0
        for offset in range(diaglen):
            if diagonal[offset] >= threshold:
                runlen += 1

            else:
                try:
                    run[runlen] += 1
                except KeyError:
                    # runlen key doesn't exist yet
                    run[runlen] = 1

                runlen = 0
                nrun += 1

        if runlen:
            try:
                run[runlen] += 1
            except KeyError:
                # runlen key doesn't exist yet
                run[runlen] = 1
            nrun += 1

        # insert into data frame, the dataframe is randomly ordered
        for r in run:
            try:
                i = runframe['len'].index(r)
                runframe['count'][i] += run[r]

            except ValueError:
                runframe['len'].append(r)
                runframe['count'].append(run[r])

        return len(runframe['len'])

    def sortFrame(self, frame, keyvar):
        """-----------------------------------------------------------------------------------------
        Sort all the variables in the dataframe according to the order of keyvar
        TODO should this and return the min and max values?

        :param frame: string
        :param keyvar: string
        :return: True
        -----------------------------------------------------------------------------------------"""
        unsorted = self.frame[frame]

        # save the order so it can be applied to all viariables in the dataframe
        order = sorted(range(len(unsorted[keyvar])), key=lambda x: unsorted[keyvar][x])

        sorted_frame = {}
        for column in unsorted:
            sorted_frame[column] = []
            for i in order:
                sorted_frame[column].append(unsorted[column][i])

        self.frame[frame] = sorted_frame
        return True

    def bdot(self, dataname, figurename, width=1, color=1, mode='dot', set_colormap=True):
        """-----------------------------------------------------------------------------------------
        Bokeh plot of dots in the main panel, and colorbar in the legend panel

        :param dataname: string, name of a dataframe in self.frame
        :param figurename: string, a figure defined in setupBokeh and stored in self.figure
        :param width: boolean, scale size of markers by the score
        :param color: boolean, scale the color of the markers by the score
        :param mode: string, if dot use the circle renderer, otherwise segment renderer
        :param set_colormap: boolean, set the colormap based on score range, turn off for second
        plot to use the same scale
        :return: True
        -----------------------------------------------------------------------------------------"""
        data = self.frame[dataname]
        figure = self.figure[figurename]
        legend = self.figure['legend']
        window = self.window
        threshold = self.threshold
        alpha = self.alpha

        scoremin, scoremax = self.valueMinMax(data['score'])

        if width == 1:
            self.scaleColumn('dots', 'score', 'size',
                             (threshold - 1, scoremax), (self.mindotsize, self.maxdotsize))
        else:
            data['size'] = [self.mindotsize for _ in range(len(data['score']))]

        if color == 1:
            pass
        else:
            data['score'] = [scoremax for _ in range(len(data['score']))]

        if set_colormap:
            if color == 1:
                cmap = LinearColorMapper(self.palette,
                                         low=max(threshold - 1.0, scoremin - 1), high=scoremax)
            else:
                cmap = LinearColorMapper(self.palette,
                                         low=threshold-0.1, high=threshold)
            self.cmap = cmap
        else:
            cmap = self.cmap

        source = ColumnDataSource(data)
        if mode == 'dot':
            # figure.circle(source=source, x='x', y='y', size='size',
            #               line_color=transform('score', cmap), line_alpha=alpha,
            #               fill_color=transform('score', cmap), fill_alpha=alpha)
            glyph = Scatter(x='x', y='y', size='size',
                          line_color=transform('score', cmap), line_alpha=alpha,
                          fill_color=transform('score', cmap), fill_alpha=alpha)
            figure.add_glyph(source, glyph)

        else:
            # line mode
            figure.segment(source=source, x0='x', x1='x1', y0='y', y1='y1', line_width='size',
                           line_color=transform('score', cmap), alpha=alpha)

        # color bar is in a separate window, self.legend, so it doesn't disturb the
        # aspect ratio
        if color:
            color_bar = ColorBar(color_mapper=cmap, label_standoff=3, bar_line_color='black',
                                 scale_alpha=alpha, width=20, margin=0, location=(0, 0),
                                 major_tick_in=20, major_tick_out=5, major_tick_line_color='black')

            legend.add_layout(color_bar, 'left')

        return True

    def bscoreDist(self, figurename, dataname, color):
        """-----------------------------------------------------------------------------------------
        Bokeh plot of score distribution and cumulative score distribution.

        :param figurename: string, name of figures (stored in self.figure)
        :param dataname: string, name of data frame (stored in self.frame)

        :param figurename: string, name of figures (stored in self.figure)
        :param dataname: string, name of data frame (stored in self.frame)
        :param color: string, and valid Bokeh color, used to fill bars
        :return: True
        -----------------------------------------------------------------------------------------"""
        data = self.frame[dataname]
        figure = self.figure[figurename]

        minp, maxp = self.valueMinMax(data['count'])

        source = ColumnDataSource(data)

        # observed score density
        figure.vbar(source=source, x='score', top='count', width=0.8, color=color,
                    line_color='black', alpha=self.alpha, bottom=0.0)
        figure.y_range = Range1d(0.0, maxp * 1.1)

        return True

    def brunDist(self, figurename, dataname, color):
        """-----------------------------------------------------------------------------------------
        Bokeh plot of run length distribution

        :param figurename: string, name of figures (stored in self.figure)
        :param dataname: string, name of data frame (stored in self.frame)
        :param color: string, and valid Bokeh color, used to fill bars
        :return: True
        -----------------------------------------------------------------------------------------"""
        run = self.frame[dataname]
        figure = self.figure[figurename]

        source = ColumnDataSource(run)
        minrun = 1
        # x = [i for i in range(minrun, maxrun + 1)]
        # observed and simulated run lengths,  need bottom=1 because of log axis

        figure.vbar(source=source, x='len', top='count', width=0.8,
                    color=color, line_color='black', alpha=self.alpha,
                    line_width=0.5, bottom=0.1)

        return True

    def bscoreCumulative(self, figurename, dataname):
        """-----------------------------------------------------------------------------------------
        Bokeh plot of cumulative distribution as a line on right hand axis

        :param figurename: string, name of figures (stored in self.figure)
        :param dataname: string, name of data frame (stored in self.frame)
        :return: True
        -----------------------------------------------------------------------------------------"""
        data = self.frame[dataname]
        figure = self.figure[figurename]

        source = ColumnDataSource(data)

        figure.extra_y_ranges = {"cumulative": Range1d(start=0.0, end=1.0)}
        axis2 = LinearAxis(y_range_name="cumulative")
        axis2.ticker.num_minor_ticks = 10
        figure.add_layout(axis2, 'right')
        figure.line(source=source, x='score', y='cumulative', y_range_name='cumulative',
                    line_width=2, color='#1122cc')

        # shaded box showing 95% level
        box = BoxAnnotation(bottom=0.95, top=1.0, y_range_name='cumulative',
                            fill_color='#FFBBBB', line_width=3, line_dash='dashed')
        figure.add_layout(box)

        return True

    def writeFrame(self, framename, key='x', out=sys.stdout):
        """-----------------------------------------------------------------------------------------
        Write the dataframe out as a table to the specified output file.  Output file should be
        opened for writing in advance.

        TODO figure out how to format values more nicely

        :param framename: string, name of a dataframe in self.frame
        :param key: string, name of column to use as key (first column in table)
        :param out: open output file
        :return: True
        -----------------------------------------------------------------------------------------"""
        frame = self.frame[framename]

        out.write('\n{} dataframe\n'.format(framename))
        out.write('\t{}'.format(key))
        for column in frame:
            if column == key:
                continue
            out.write('\t{}'.format(column))
        out.write('\n')

        n = len(frame[key])
        for i in range(n):
            out.write('\t{}'.format(frame[key][i]))
            for column in frame:
                if column == key:
                    continue
                out.write('\t{}'.format(frame[column][i]))
            out.write('\n')

        return True

    def show(self, *args, **kwargs):
        """-----------------------------------------------------------------------------------------
        Delegate to plt.show().  Makes syntax a little easier in application since the object is
        used instead of the plotting class

        :param args: arguments to pass to show()
        :param kwargs: arguments to pass to show()
        :return: True
        -----------------------------------------------------------------------------------------"""
        show(self.grid, *args, **kwargs)

        return True

    @staticmethod
    def cumulative(score, total):
        """-----------------------------------------------------------------------------------------
        Return cumulative score probability distribution as a list.

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

    def addCumulative(self, data, sourcecol, destcol):
        """-----------------------------------------------------------------------------------------
        Add cumulative distribution to dataframe data, based on column sourcecol and stored in a
        new column named destcol

        :param data: string (dataframe in self.frames)
        :param sourcecol: string, column name in self.frame[data]
        :param destcol: string, new column name for cumulative distribution
        :return: True
        -----------------------------------------------------------------------------------------"""
        data = self.frame[data]
        source = data[sourcecol]
        cum = []
        total = 0
        for v in source:
            total += v
            cum.append(total)

        for i in range(len(cum)):
            cum[i] /= total

        data[destcol] = cum
        return True

    def addSegment(self, framename, xcol='x', ycol='y', xnew='x1', ynew='y1'):
        """-----------------------------------------------------------------------------------------
        convert x, y dot positions to line segments; the segment renderer requires beginning and
        ending points for each segment. The existing x and y are modified to be the beginning and
        new variables (xnew and ynew) are added for the end points.
        
        :param xcol: string, name of x column in data frame
        :param ycol: string, name of y column in data frame
        :param xnew: string, name of new x column in data frame (end of segment)
        :param ynew: string, name of new y column in data frame (end of segment)
        :return: True
        -----------------------------------------------------------------------------------------"""
        frame = self.frame[framename]
        frame[xnew] = []
        frame[ynew] = []

        # correct the direction when sequence 2 is reversed
        yinc = self.yinc
        dither = 0.5
        ydither = [-dither * yinc, dither * yinc]

        for pos in range(len(frame[xcol])):
            frame[xnew].append(frame[xcol][pos] + dither)
            frame[ynew].append(frame[ycol][pos] + ydither[1])
            frame[xcol][pos] -= dither
            frame[ycol][pos] += ydither[0]

        return True

    @staticmethod
    def density(score, total):
        """-----------------------------------------------------------------------------------------
        Convert a list representing the score distribution to a density by dividing by total
        :param score: list of int or float
        :param total: total number of scores (sum(score))
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

    @staticmethod
    def valueMinMax(score):
        """-----------------------------------------------------------------------------------------
        Returns the minimum and maximum value in a list of values.

        :param score: list
        :return: float, float
        -----------------------------------------------------------------------------------------"""
        if score:
            scoremin = score[0]
            scoremax = score[0]
        else:
            scoremin = scoremax = 0

        for s in score:
            scoremin = min(scoremin, s)
            scoremax = max(scoremax, s)

        return scoremin, scoremax


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    match = Diagonal()

    fasta = Fasta(filename=sys.argv[1])
    fasta.read()
    fasta1 = fasta.copy()

    fasta2 = fasta.copy()
    fasta2.id = 'seq 2'
    fasta2.doc = 'Sequence 2'

    fasta1.seq = fasta1.seq[:200]
    fasta2.seq = fasta2.seq[:400]

    dataframes = [{'data': 'dots', 'fn': match.windowThreshold, 'var': ['x', 'y', 'score']},
                  {'data': 'scoredist', 'fn': match.histogramScore, 'var': ['score', 'count']},
                  {'data': 'rundist', 'fn': match.histogramRun, 'var': ['len', 'count']},
                  {'data': 'randomscore', 'fn': None, 'var': ['score', 'count']},
                  {'data': 'randomrun', 'fn': None, 'var': ['len', 'count']}
                  ]

    # select list of tests to run, one plot in browser for each
    # tests = [0,1,2,3]
    tests = range(8)
    # tests = [2]

    for test in tests:

        if test == 0:
            w = 4
            t = 1
            match.title = 'test {} - forward dotplot size cueing only, w={} t={}'. \
                format(test, w, t)
            match.setupFrame(dataframes)
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=True, color=False)

        elif test == 1:
            w = 4
            t = 1
            match.title = 'test {} - forward dotplot color cueing only, w={} t={}'. \
                format(test, w, t)
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=False, color=True)

        elif test == 2:
            w = 12
            t = 8
            match.title = 'test {} - forward dotplot size, color, stats, w={} t={}'. \
                format(test, w, t)

            fasta1 = fasta.copy()
            fasta2 = fasta.copy()

            match.setupFrame(dataframes)
            match.setupCalculation(fasta1, fasta1, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=True, color=True)

            match.sortFrame('scoredist', 'score')
            match.addCumulative('scoredist', 'count', 'cumulative')
            match.bscoreDist('scoredist', 'scoredist', color='#0000ff')
            match.bscoreCumulative('scoredist', 'scoredist')
            match.brunDist('rundist', 'rundist', color='#0000ff')

        elif test == 3:
            w = 14
            t = 9
            match.title = 'test {} - forward/reverse dotplot, w={} t={}'. \
                format(test, w, t)

            fasta1 = fasta.copy()
            fasta2 = fasta.copy()

            match.setupFrame(dataframes)
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=True, color=True)

            match.seqreverse = True
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t, resetstat=False)
            # match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=True, color=True, set_colormap=False)

            match.sortFrame('scoredist', 'score')
            match.addCumulative('scoredist', 'count', 'cumulative')
            match.bscoreDist('scoredist', 'scoredist', color='#0000ff')
            match.bscoreCumulative('scoredist', 'scoredist')
            match.brunDist('rundist', 'rundist', color='#0000ff')

            match.seqreverse = False  # for next plot
            match.yinc = 1

        elif test == 4:
            w = 16
            t = 9
            match.title = 'test {} - forward/reverse dotplot with random , w={} t={}'. \
                format(test, w, t)

            fasta1 = fasta.copy()
            fasta2 = fasta.copy()

            match.setupFrame(dataframes)
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=True, color=True)

            match.seqreverse = True
            match.resetFrame('dots')
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t, resetstat=False)
            # match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=True, color=True, set_colormap=False)

            match.sortFrame('scoredist', 'score')
            match.addCumulative('scoredist', 'count', 'cumulative')
            match.bscoreDist('scoredist', 'scoredist', color='#0000ff')
            match.bscoreCumulative('scoredist', 'scoredist')
            match.brunDist('rundist', 'rundist', color='#0000ff')

            # random score distribution
            match.single = True
            match.random(n=match.nscore)
            match.histogramScore('randomscore', 1)
            match.sortFrame('randomscore', 'score')
            match.bscoreDist('scoredist', 'randomscore', color='#ff0000')

            match.histogramRun('randomrun', 1)
            match.brunDist('rundist', 'randomrun', color='#ff0000')
            match.single = False

            match.writeFrame('scoredist', key='score')
            match.writeFrame('rundist', key='len')

            match.seqreverse = False  # for next plot
            match.yinc = 1

        elif test == 5:
            w = 16
            t = 9
            match.title = 'test {} - forward/reverse lineplot with random , w={} t={}'. \
                format(test, w, t)

            fasta1 = fasta.copy()
            fasta2 = fasta.copy()

            match.setupFrame(dataframes)
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.addSegment('dots')
            match.bdot('dots', 'main', width=True, color=True, mode='line')

            match.seqreverse = True
            match.resetFrame('dots')
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t, resetstat=False)
            # match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.addSegment('dots')
            match.bdot('dots', 'main', width=True, color=True, mode='line', set_colormap=False)

            match.sortFrame('scoredist', 'score')
            match.addCumulative('scoredist', 'count', 'cumulative')
            match.bscoreDist('scoredist', 'scoredist', color='#0000ff')
            match.bscoreCumulative('scoredist', 'scoredist')
            match.brunDist('rundist', 'rundist', color='#0000ff')

            # random score distribution
            match.single = True
            match.random(n=match.nscore)
            match.histogramScore('randomscore', 1)
            match.sortFrame('randomscore', 'score')
            match.bscoreDist('scoredist', 'randomscore', color='#ff0000')

            match.histogramRun('randomrun', 1)
            match.brunDist('rundist', 'randomrun', color='#ff0000')
            match.single = False

            match.writeFrame('scoredist', key='score')
            match.writeFrame('rundist', key='len')

            match.seqreverse = False  # for next plot
            match.yinc = 1

        elif test == 6:
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

            match.setupFrame(dataframes)
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.bdot('dots', 'main', width=True, color=True)

            match.sortFrame('scoredist', 'score')
            match.addCumulative('scoredist', 'count', 'cumulative')
            match.bscoreDist('scoredist', 'scoredist', color='#0000ff')
            match.bscoreCumulative('scoredist', 'scoredist')
            match.brunDist('rundist', 'rundist', color='#0000ff')

            # random score distribution
            match.single = True
            match.random(n=match.nscore)
            match.histogramScore('randomscore', 1)
            match.sortFrame('randomscore', 'score')
            match.bscoreDist('scoredist', 'randomscore', color='#ff0000')

            match.histogramRun('randomrun', 1)
            match.brunDist('rundist', 'randomrun', color='#ff0000')
            match.single = False

            match.writeFrame('scoredist', key='score')
            match.writeFrame('rundist', key='len')

        elif test == 7:
            w = 16
            t = 10
            match.title = 'test {} - protein lineplot with blosum, w={} t={}'. \
                format(test, w, t)
            match.readNCBI('table/BLOSUM62.matrix')

            fasta = Fasta(filename='../data/globin.pfa')
            fasta.read()
            fasta1 = fasta.copy()
            fasta.read()
            fasta2 = fasta.copy()

            match.setupFrame(dataframes)
            match.setupCalculation(fasta1, fasta2, window=w, threshold=t)
            match.setupBokeh(cbase='Viridis', clevels=256, creverse='True')
            # match.resetFrame('scoredist')
            match.allDiagonals(select=['dots', 'scoredist', 'rundist'])
            match.addSegment('dots')
            match.bdot('dots', 'main', width=True, color=True, mode='line')

            match.sortFrame('scoredist', 'score')
            match.addCumulative('scoredist', 'count', 'cumulative')
            match.bscoreDist('scoredist', 'scoredist', color='#0000ff')
            match.bscoreCumulative('scoredist', 'scoredist')
            match.brunDist('rundist', 'rundist', color='#0000ff')

            # random score distribution
            match.single = True
            match.random(n=match.nscore)
            match.histogramScore('randomscore', 1)
            match.sortFrame('randomscore', 'score')
            match.bscoreDist('scoredist', 'randomscore', color='#ff0000')

            match.histogramRun('randomrun', 1)
            match.brunDist('rundist', 'randomrun', color='#ff0000')
            match.single = False

            match.writeFrame('scoredist', key='score')
            match.writeFrame('rundist', key='len')

        match.show()

    exit(0)

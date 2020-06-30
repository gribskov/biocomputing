"""=================================================================================================
Even though RLE provides a sparse matrix representation, it is difficult to save scores for
positions in the matching windows.  As an alternative calculate and plot one diagonal at a time

TODO: plot score distribution
TODO: plot run length distribution (filtered)
TODO: add correct units to cmap display


Michael Gribskov     25 June 2020
================================================================================================="""
import sys
from datetime import date
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
        self.fig = plt.figure(figsize=[11.0, 8.5])
        self.ax = None
        self.title = ''
        self.threshold = 0
        self.window = 0

        self.s1 = Fasta()
        self.s2 = Fasta()
        self.i1 = None  # integer array representation of sequences
        self.i2 = None
        self.l1 = 0
        self.l2 = 0

    def setup(self, seq1, seq2, window=5, threshold=3):
        """-----------------------------------------------------------------------------------------
        Load the sequences and do some basic setup for score calculations and plotting. Sequences
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

        # plot setup
        if self.title:
            titlestr = self.title
        else:
            now = date.today()
            titlestr = 'Dotplot of {} and {} - {}'.format(self.s1.id, self.s2.id, now)

        self.fig.suptitle(titlestr)
        if not self.ax:
            self.ax = self.fig.add_subplot(1, 1, 1)

        ax = self.ax
        ax.set_aspect(1.0)

        ax.set_xlim(0, self.l1 + 1)
        ax.set_ylim(0, self.l2 + 1)

        ax.set_xlabel('\n'.join([self.s1.id, self.s1.doc]))
        ax.set_ylabel('\n'.join([self.s2.doc, self.s2.id]))

        return True

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

    # def drawDot(self, rev=False, width=False, color=False):
    #     """-----------------------------------------------------------------------------------------
    #     Draw all diagonals as dots.  the score is reported in the first position of the window so
    #     the position of the dot must be offset to lie in the middle of the window.  Zero origin
    #     coordinates must also be incremented by 1
    # 
    #     :return: True
    #     -----------------------------------------------------------------------------------------"""
    #     window = self.window
    #     halfwindow = (window - 1) / 2.0
    #     threshold = self.threshold
    #     diagonal = self.diagonal
    #     symbol = 'ko'
    #     l2 = self.l2
    # 
    #     # markersize scaled to score
    #     minmarker = 0.1
    #     maxmarker = 6.0
    #     defmarker = 2.0
    #     woff = 0
    #     wscale = 1
    #     if width:
    #         wscale = (maxmarker - minmarker) / (window - threshold)
    #         woff = 1 - wscale * threshold
    # 
    #     # reversed plot: invert the direction and change color
    #     yinc = 1
    #     if rev:
    #         yinc = -1
    #         symbol = 'r.'
    # 
    #     for d in range(self.l1 + self.l2 - 1):
    #         dscore = self.diagonalScore(d)
    #         diaglen, xpos, ypos = self.diagLenBegin(d)
    #         xpos += halfwindow
    #         if rev:
    #             ypos = l2 - ypos - halfwindow - 1
    #         else:
    #             ypos += halfwindow
    # 
    #         for pos in range(diaglen - window + 1):
    #             size = defmarker
    #             if dscore[pos] >= threshold:
    #                 if width:
    #                     size = dscore[pos] * wscale + woff
    # 
    #                 plt.plot(xpos, ypos, symbol, markersize=size)
    # 
    #             xpos += 1
    #             ypos += yinc
    # 
    #     return True

    def drawDot(self, rev=False, cmap='gray', colreverse=True, width=False, color=False):
        """-----------------------------------------------------------------------------------------
        Draw all diagonals as dots.  the score is reported in the first position of the window so
        the position of the dot must be offset to lie in the middle of the window.  Zero origin
        coordinates must also be incremented by 1

        :return: True
        -----------------------------------------------------------------------------------------"""
        window = self.window
        halfwindow = (window - 1) / 2.0
        threshold = self.threshold
        msize = 4.0
        diagonal = self.diagonal
        l2 = self.l2

        cc = plt.cm.get_cmap(cmap, window - threshold)
        if colreverse:
            cc = cc.reversed()
        plt.colorbar(cm.ScalarMappable(cmap=cc), shrink=0.4, fraction=0.05)
        self.ax.set_facecolor('w')
        # self.ax.set_facecolor(cc(0.0))

        # reversed plot: invert the direction and change color
        # markers are slightly smaller on reversed plot
        yinc = 1
        if rev:
            yinc = -1
            msize /= 1.5
            # symbol = 'r.'

        maxmarker = 1
        minmarker = 0
        wscale = 1.0 / (window - threshold)
        woff = (1 / (window - threshold) - wscale * threshold)

        for d in range(self.l1 + self.l2 - 1):
            dscore = self.diagonalScore(d)
            diaglen, xpos, ypos = self.diagLenBegin(d)
            xpos += halfwindow
            if rev:
                ypos = l2 - ypos - halfwindow - 1
            else:
                ypos += halfwindow

            for pos in range(diaglen - window + 1):
                if dscore[pos] >= threshold:
                    # f=str(1.0-dscore[pos]/window)
                    f = dscore[pos] * wscale + woff
                    size = msize
                    if width:
                        size *= f
                    if not color:
                        f = 1.0
                    plt.plot(xpos, ypos, 'o', markersize=size, c=cc(f), alpha=0.75)

                xpos += 1
                ypos += yinc

        return True

    def drawLine(self, dither=0.05):
        """-----------------------------------------------------------------------------------------
        Draw all diagonals as lines.  the score is reported in the first position of the window so
        the position of the dot must be offset to lie in the middle of the window.  Zero origin
        coordinates must also be incremented by 1

        :return: True
        -----------------------------------------------------------------------------------------"""
        window = self.window
        threshold = self.threshold
        diagonal = self.diagonal
        halfwindow = window / 2.0 + 1
        dither2 = dither * 2

        for d in range(self.l1 + self.l2 - 1):
            dscore = self.diagonalScore(d)
            diaglen, xpos, ypos = self.diagLenBegin(d)
            xpos += halfwindow + dither
            ypos += halfwindow + dither

            line = 0
            begin = []
            for pos in range(diaglen - window + 1):
                if dscore[pos] >= threshold:
                    # this window is on a line
                    # if already on a line, do nothing until end of line
                    if not line:
                        # not in a line, start a new line
                        begin[0:1] = [xpos - dither2, ypos - dither2]

                    line += 1

                else:
                    # this window is not part of line, draw the line
                    if line:
                        plt.plot([begin[0], xpos - 1], [begin[1], ypos - 1], color='k')
                        line = 0

                xpos += 1
                ypos += 1

            if line:
                # ended in a line, draw the last one
                plt.plot([begin[0], xpos - 1], [begin[1], ypos - 1], color='k')

        return True

    def drawLineWidth(self):
        """-----------------------------------------------------------------------------------------
        Draw all diagonals as dots.  the score is reported in the first position of the window so
        the position of the dot must be offset to lie in the middle of the window.  Zero origin
        coordinates must also be incremented by 1

        :return: True
        -----------------------------------------------------------------------------------------"""

        minmarker = 0.5
        maxmarker = 6.0

        window = self.window
        halfwindow = window / 2.0 + 1
        threshold = self.threshold
        diagonal = self.diagonal

        for d in range(self.l1 + self.l2 - 1):
            dscore = self.diagonalScore(d)
            diaglen, xpos, ypos = self.diagLenBegin(d)
            xpos += halfwindow
            ypos += halfwindow
            for pos in range(diaglen - window + 1):
                if dscore[pos] >= threshold:
                    size = (dscore[pos] - threshold + 1) / (window - threshold + 1) * (
                            maxmarker - minmarker)
                    # f = (dscore[pos] - threshold) / (window - threshold)
                    plt.plot([xpos - 0.5, xpos + 0.5], [ypos - 0.5, ypos + 0.5], color='k', lw=size,
                             solid_capstyle='round')

                xpos += 1
                ypos += 1

        return True

    def drawLineColor(self, cmap='gray', reverse=True):
        """-----------------------------------------------------------------------------------------
        Draw all diagonals as dots.  the score is reported in the first position of the window so
        the position of the dot must be offset to lie in the middle of the window.  Zero origin
        coordinates must also be incremented by 1

        :return: True
        -----------------------------------------------------------------------------------------"""
        window = self.window
        halfwindow = window / 2.0 + 1
        threshold = self.threshold
        diagonal = self.diagonal

        cc = plt.cm.get_cmap(cmap, window + 1 - threshold)
        if reverse:
            cc = cc.reversed()
        plt.colorbar(cm.ScalarMappable(cmap=cc), fraction=0.1, shrink=0.4, values=range(
            threshold, window + 1))
        self.ax.set_facecolor(cc(0.0))

        for d in range(self.l1 + self.l2 - 1):
            dscore = self.diagonalScore(d)
            diaglen, xpos, ypos = self.diagLenBegin(d)
            xpos += halfwindow
            ypos += halfwindow
            for pos in range(diaglen - window + 1):
                if dscore[pos] >= threshold:
                    # f=str(1.0-dscore[pos]/window)
                    f = (dscore[pos] - threshold) / (window - threshold)
                    plt.plot([xpos - 0.5, xpos + 0.5], [ypos - 0.5, ypos + 0.5], c=cc(f), lw=2,
                             solid_capstyle='round')

                xpos += 1
                ypos += 1

        return True

    def show(self, *args, **kwargs):
        """-----------------------------------------------------------------------------------------
        Delegate to plt.show().  Makes syntax a little easier in application since the object is
        used instead of the matplotlib class (application doesn't have to know matplotlib).

        :param args:
        :param kwargs:
        :return: True
        -----------------------------------------------------------------------------------------"""
        plt.show(*args, **kwargs)


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

    # match.setup(fasta1, fasta2)
    match.setup(fasta1, fasta2, window=20, threshold=8)
    # match.drawDot(width=True)
    # fasta2.seq = fasta2.reverseComplement()
    # match.setup(fasta1, fasta2, window=10, threshold=6)
    # match.drawDot(rev=True, width=True)
    # # match.drawDot(rev=True, width=True)
    match.drawDot(cmap='hot', colreverse=True, color=True)

    # match.drawLineWidth()
    # match.drawDot(cmap='viridis', colreverse=False)
    # match.drawDot(cmap='Blues', colreverse=False, color=True, width=True)
    # fasta2.seq = fasta2.reverseComplement()
    # match.setup(fasta1, fasta2, window=9, threshold=6)
    # match.drawDot(cmap='Reds', rev=True, colreverse=False, width=True)

    # match.drawLine()
    # match.drawLineColor(cmap='hot', colreverse=True)
    match.show()

    exit(0)

"""=================================================================================================
TODO need to offset lines by 1/2 window

Michael Gribskov     24 June 2020
================================================================================================="""
from sequence.score import Score
from wordmatch import Base


class Windowmatch(Base, Score):
    """=============================================================================================
    Base provides some calculations related to diagonals and Fasta sequences
    Score provides scoring table methods

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Constructor for class windowmatch.  Sequences s1 and s2 are defined in Base class.

        -----------------------------------------------------------------------------------------"""
        Base.__init__(self)
        Score.__init__(self)
        self.diagonal = []
        self.window = 0
        self.threshold = 0

    def seqToInt(self):
        """-----------------------------------------------------------------------------------------
        Convert sequence strings to an integer arrays and stores in object.  An integer array is
        more convenient for direct lookups in the scoring table than a string

        :return: int, int length of sequence listss
        -----------------------------------------------------------------------------------------"""
        a2i = self.a2i

        self.i1 = [a2i[c] for c in self.s1.seq]
        self.i2 = [a2i[c] for c in self.s2.seq]

        return len(self.i1), len(self.i2)

    def windowScore(self):
        """-----------------------------------------------------------------------------------------
        Calculate the moving window sum of comparison score along diagonals and stores in the
        object.
        RLE runs stored in diagonal as [offset, length] where offset is the end position of the
        windown along the diagonal (same as wordmatch).

        :return:
        -----------------------------------------------------------------------------------------"""
        l1, l2 = self.seqToInt()
        i1 = self.i1
        i2 = self.i2

        window = self.window
        threshold = self.threshold
        cmp = self.table

        diagonal = self.diagonal
        nmatch = 0
        for d in range(l1 + l2 - 1):
            diagonal.append([])
            diaglen, pos1, pos2 = Base.diagLenBegin(d, l1, l2)
            old1 = pos1
            old2 = pos2

            if diaglen < window:
                # skip   diagonals shorter than window length
                continue

            runlen = 0
            score = 0

            # first window
            for offset in range(window):
                score += cmp[i1[pos1]][i2[pos2]]
                pos1 += 1
                pos2 += 1

            if score >= threshold:
                runlen = 1
                nmatch += 1

            # rest of diagonal
            for offset in range(window, diaglen):
                score -= cmp[i1[old1]][i2[old2]]
                score += cmp[i1[pos1]][i2[pos2]]

                if score >= threshold:
                    runlen += 1
                    nmatch += 1
                    nmatch += 1

                else:
                    if runlen:
                        # run ends at previous position
                        self.diagonal[d].append([offset - 1, runlen])
                        runlen = 0

                old1 += 1
                old2 += 1
                pos1 += 1
                pos2 += 1

            if runlen:
                # in a run at end of diagonal (do not subtract because the end position is
                # the true end)
                self.diagonal[d].append([offset, runlen])

        return nmatch


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    from sequence.fasta import Fasta
    from plotter import Plotter

    match = Windowmatch()
    print('done {}'.format(type(match)))
    print(match.alphabet)

    # match.readNCBI('table/NUC4.4.matrix')
    print(match.format())

    fasta1 = Fasta(filename=sys.argv[1])
    fasta1.read()

    fasta2 = Fasta()
    fasta2.id = 'seq2'
    fasta2.doc = ' bases 1:50'
    fasta2.seq = fasta1.seq[:50]

    fasta1.seq = fasta1.seq[:200]

    match.s1 = fasta1
    match.s2 = fasta2
    l1, l2 = match.seqToInt()
    print(l1, l2)

    match.window = 10
    match.threshold = 5
    nmatch = match.windowScore()
    print('window: {}     threshold: {}     nmatch: {}'. \
          format(match.window, match.threshold, nmatch))
    plot = Plotter()
    plot.match = match
    plot.setup()
    plot.lines(dither=0.03)
    # plot.dots()
    plot.show()

    exit(0)

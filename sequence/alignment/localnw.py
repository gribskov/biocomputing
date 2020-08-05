"""=================================================================================================
just a trial of using an object to represent the cell and pointers.  since it is a local alignment
only non'zero cells need to be represented.  maybe this will make it sparse enough to be worthwhile

Michael Gribskov     04 August 2020
================================================================================================="""
import sys
from math import log10
import random
from scipy import stats
from sequence.fasta import Fasta
from sequence.score import Score


class Cell:
    """=============================================================================================

    ============================================================================================="""
    count = 0

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.n = Cell.count
        Cell.count += 1
        self.score = 0
        self.p = []


class Alignment(Score):
    """---------------------------------------------------------------------------------------------

    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        Score.__init__(self)
        self.s1 = None
        self.s2 = None
        self.i1 = None
        self.i2 = None

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

    def localBrute(self, open, extend):
        """-----------------------------------------------------------------------------------------
        Local alignment score only.
        s1 is the horizontal sequence and s2 is the vertical sequence.  this makes s2 the row
        index and s1 the column index.

        :param open: float, gap opening penalty
        :param extend: float, gap extension penalty
        :return:
        -----------------------------------------------------------------------------------------"""
        cmp = self.table
        i1 = self.i1
        i2 = self.i2
        l1 = len(i1)
        l2 = len(i2)

        edge = Cell()  # a dummy cell for the edges
        edge.score = self.min + open

        scoremax = 0
        posmax = [0, 0]
        score = [[Cell() for i in i1] for j in i2]
        bestrow = Cell()
        bestrow.p = edge
        bestrow.score = edge.score
        bestcol = []
        for _ in i1:
            c = Cell()
            c.p = edge
            c.score = edge.score
            bestcol.append(c)

        diag = Cell()

        jpos = 0
        for j in i2:
            # print(score) # uncomment to print score matrix
            diag.p = edge
            diag.score = edge.score
            bestrow.p = edge
            bestrow.score = edge.score + extend
            ipos = 0
            for i in i1:
                previous = max(diag.score, bestcol[ipos].score, bestrow.score, 0.0)
                cell = max(cmp[j][i] + previous, 0.0)
                if cell > 0:
                    score[jpos][ipos].score = cell
                    if cell > scoremax:
                        scoremax = cell
                        posmax = [ipos, jpos]

                    if previous > 0:
                        # only set pointers if cell score > zero and there is a non-zero best
                        # previous score
                        for dir in (diag, bestrow, bestcol[ipos]):
                            if dir.score == previous:
                                # set pointers for all directions
                                score[jpos][ipos].p.append( dir.p )

                    if diag.score + open > bestrow.score + extend:
                        # what if scores are equal? some paths missed
                        bestrow.p = diag.p
                        bestrow.score = diag.score + open

                    if diag.score + open > bestcol[ipos].score + extend:
                        bestcol[ipos].p = diag.p
                        bestcol[ipos].score = diag.score + open

                diag.p = score[jpos][ipos]
                diag.score = score[jpos][ipos].score
                score[jpos][ipos].score = cell
                ipos += 1

            jpos += 1

        # print(score)
        return scoremax, posmax


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    align = Alignment()
    #
    # align.s1 = Fasta(filename=sys.argv[1])
    # align.s2 = Fasta(filename=sys.argv[2])
    # align.readNCBI('../../dotplot/table/BLOSUM62.matrix')

    # testing
    align.s1 = Fasta()
    align.s1.seq = 'ACTGCC'
    align.s2 = Fasta()
    align.s2.seq = 'ATGCC'
    align.readNCBI('../../dotplot/table/NUC4.4.matrix')

    align.seqToInt()
    # # random.shuffle(align.i1)          # uncomment to test scores for random alignments
    original_score, bestpos = align.localBrute(-10, -1)

    # print('original score: {} at {}'.format(original_score, bestpos))

    exit(0)

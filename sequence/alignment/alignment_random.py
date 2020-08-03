"""=================================================================================================
Significance of alignment using random simulation

Michael Gribskov     02 August 2020
================================================================================================="""
import sys
from sequence.fasta import Fasta
from sequence.score import Score


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

    def localScore(self, open, extend):
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

        edge = self.min + open  # a small value to use for gaps on the edges

        bestrow = edge
        bestcol = [edge for _ in i1]
        score = [0 for _ in i1]
        cell = 0
        diag = 0

        # initialize row zero
        j = i2[0]
        pos = 0
        for i in i1:
            s = max(cmp[j][i], 0.0)
            score[pos] = s
            bestcol[pos] += extend
            pos += 1

        for j in i2[1:]:
            print(score)
            i = i1[0]
            score[0] = max(cmp[j][i], 0.0)

            diag = score[0]
            bestcol[0] = max(diag + open, bestcol[0] + extend)
            bestrow = edge + extend
            pos = 1
            for i in i1[1:]:
                cell = max(cmp[j][i] + max(diag, bestcol[pos], bestrow), 0.0)
                bestrow = max(diag + open, bestrow + extend)
                bestcol[pos] = max(diag + open, bestcol[pos] + extend)
                diag = score[pos]
                score[pos] = cell
                pos += 1

        print(score)
        return True


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    align = Alignment()

    # align.s1 = Fasta(filename=sys.argv[1])
    # align.s2 = Fasta(filename=sys.argv[2])
    # align.readNCBI('../../dotplot/table/BLOSUM62.matrix')

    # testing
    align.s1 = Fasta()
    align.s1.seq = 'AATGCC'
    align.s2 = Fasta()
    align.s2.seq = 'AAGCC'
    align.readNCBI('../../dotplot/table/NUC4.4.matrix')

    align.seqToInt()
    align.localScore(-1, -1)

    exit(0)

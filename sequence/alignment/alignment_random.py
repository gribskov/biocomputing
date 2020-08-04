"""=================================================================================================
Significance of alignment using random simulation

Michael Gribskov     02 August 2020
================================================================================================="""
import sys
from math import log10
import random
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

        scoremax = 0
        bestrow = edge
        bestcol = [edge for _ in i1]
        score = [0 for _ in i1]
        cell = 0
        diag = 0

        for j in i2:
            # print(score)
            diag = 0  # uncomment to print score matrix
            bestrow = edge + extend
            pos = 0
            for i in i1:
                cell = max(cmp[j][i] + max(diag, bestcol[pos], bestrow), 0.0)
                scoremax = max(cell, scoremax)
                bestrow = max(diag + open, bestrow + extend)
                bestcol[pos] = max(diag + open, bestcol[pos] + extend)
                diag = score[pos]
                score[pos] = cell
                pos += 1

        # print(score)
        return scoremax


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    align = Alignment()

    align.s1 = Fasta(filename=sys.argv[1])
    align.s2 = Fasta(filename=sys.argv[2])
    align.readNCBI('../../dotplot/table/BLOSUM62.matrix')

    # testing
    # align.s1 = Fasta()
    # align.s1.seq = 'ACTGCC'
    # align.s2 = Fasta()
    # align.s2.seq = 'ATGCC'
    # align.readNCBI('../../dotplot/table/NUC4.4.matrix')

    align.seqToInt()
    score = align.localScore(-10, -1)
    print('original score: {}'.format(score))

    nrandom = 2000
    rscore = []
    for i in range(nrandom):
        random.shuffle(align.i2)
        # random.shuffle(align.i1)
        score = align.localScore(-20, -2)
        rscore.append(score)
        print('\t{}: {}'.format(i, score))

    i = 0
    N = []
    rscore.sort(reverse=True)
    ro = rscore[0]
    for r in rscore:
        if r < ro:
            N.append([i, ro])
            ro = r

        i += 1
        print('\t{}\t{}\t{:.3f}'.format(i, log10(i), r))
    N.append([i,ro])


    for pair in N:
        print('\t{}\t{:.3f}\t{}'.format(pair[0], log10(pair[0]), pair[1]))

    exit(0)

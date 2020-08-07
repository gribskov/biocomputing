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
        self.score = None

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
        # i dimension of score is length + 1 to avoid update exception at end of row fill loop
        score = [[Cell() for i in i1] for j in i2]
        self.score = score
        bestrow = Cell()
        # bestrow.p = edge
        # bestrow.score = edge.score
        bestcol = []
        for _ in i1:
            c = Cell()
            c.p = [edge]
            c.score = edge.score
            bestcol.append(c)

        for ipos in range(l1):
            score[0][ipos].score = edge.score
            score[0][ipos].p = [edge]

        diag = Cell()

        jpos = 0
        for j in i2:
            bestrow.p = [edge]
            bestrow.score = edge.score + extend
            diag.p = [edge]
            diag.score = edge.score

            ipos = 0
            for i in i1:
                previous = max(diag.score, bestcol[ipos].score, bestrow.score, 0.0)
                cell = max(cmp[j][i] + previous, 0.0)
                if cell > 0:
                    # score[jpos][ipos].score = cell
                    if cell > scoremax:
                        scoremax = cell
                        posmax = [jpos, ipos]

                    if previous > 0:
                        # only set pointers if cell score > zero and there is a non-zero best
                        # previous score
                        for dir in (diag, bestrow, bestcol[ipos]):
                            if dir.score == previous:
                                # set pointers for all directions
                                score[jpos][ipos].p.append(dir.p)

                # update best row and column values
                if diag.score + open > bestrow.score + extend:
                    # what if scores are equal? some paths missed
                    bestrow.p = diag.p
                    bestrow.score = diag.score + open
                else:
                    bestrow.score += extend

                if diag.score + open > bestcol[ipos].score + extend:
                    bestcol[ipos].p = diag.p
                    bestcol[ipos].score = diag.score + open
                else:
                    bestcol[ipos].score += extend

                # diagonal score for next cell
                if jpos > 0:
                    diag.p = score[jpos - 1][ipos]
                    diag.score = score[jpos - 1][ipos].score

                score[jpos][ipos].score = cell
                ipos += 1

            jpos += 1

        return scoremax, posmax

    def trace1(self, pos):
        """-----------------------------------------------------------------------------------------
        Trace back one alignment using the first pointer for each cell

        :param pos: list of 2 int, traceback start position
        :return:
        -----------------------------------------------------------------------------------------"""
        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)
        score = self.score
        current = score[pos[0]][pos[1]]

        a1 = ''
        a2 = ''
        rowold = pos[0]
        colold = pos[1]
        while current.score > 0:
            n = current.n - 1
            row = n // l1
            col = n % l1
            for c in range(col,colold-1):
                a1 += s1[c]
                a2 += '.'
            a1 += s1[col]

            for r in range(row, rowold-1):
                a1 += '.'
                a2 += s2[r]
            a2 += s2[row]

            if len(current.p):
                current = current.p[0]
            rowold = row
            colold = col

        return a1[::-1], a2[::-1]

    def writeScoreMatrix(self, file, decimal=0):
        """-----------------------------------------------------------------------------------------
        Write out the score matrix in an aligned table. Could provide scoremax, but then it wouldn't
        work for a sub-table.

        :param file: open filehandle for output, think stdout
        :param decimal:int, number of digits past decimal point
        :return: float, maximum score
        -----------------------------------------------------------------------------------------"""
        # first find the largest value (maximum column width)
        score = self.score
        scoremax = 0.0
        for row in score:
            for col in row:
                scoremax = max (scoremax, col.score)

        fmt = '{{:>{}.{}f}}'.format(len(str(scoremax))+1,decimal)

        for row in score:
            for col in row:
                file.write(fmt.format(col.score))

            file.write('\n')

        return scoremax


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
    align.s1.seq = 'ACTGCCTTGATC'
    align.s2 = Fasta()
    align.s2.seq = 'ATGCCAAAGATC'
    align.readNCBI('../../dotplot/table/NUC4.4.matrix')

    align.seqToInt()
    # # random.shuffle(align.i1)          # uncomment to test scores for random alignments
    original_score, bestpos = align.localBrute(-1, -1)
    print('original score: {} at {}'.format(original_score, bestpos))
    align.writeScoreMatrix(sys.stdout)
    a1, a2 = align.trace1(bestpos)
    print('\n{}\n{}'.format(a1,a2))

    exit(0)

from numpy.f2py.crackfortran import previous_context

from sequence.fasta import Fasta
from sequence.score import Score


class Alignment:
    """=============================================================================================
    Sequence alignment class

    ============================================================================================="""

    def __init__(self):
        self.s1 = Fasta()
        self.s2 = Fasta()
        self.i1 = []
        self.i2 = []
        self.score = Score()
        self.smat = None
        self.pmat = None
        self.gi = 0  # length independent gap penalty
        self.gd = 0  # length dependent gap penalty

    def index(self, i1=True, i2=True):
        """-----------------------------------------------------------------------------------------
        create an integer array version of the sequences.  this makes it easier to lookup scores.
        i1  and i2 indicate whether sequence 1 and sequence 2 should be indexed.

        :param self:
        :return:
        -----------------------------------------------------------------------------------------"""
        score = self.score
        if i1:
            self.i1 = []
            for letter in self.s1.seq:
                self.i1.append(score.a2i[letter])
        if i2:
            self.i2 = []
            for letter in self.s2.seq:
                self.i2.append(score.a2i[letter])

        return len(self.i1), len(self.i2)

    def matString(align, width=2):
        """-----------------------------------------------------------------------------------------
        return a string with the formatted score and path matrices
        :param width: int, field width of columns
        :return: string
        -----------------------------------------------------------------------------------------"""
        # pad the sequences with a gap at the ends
        p1 = '.' + align.s1.seq + '.'
        p2 = '.' + align.s2.seq + '.'

        s = 'score\n\n  '

        for j in range(0, len(p2)):
            s += "{} ".format(p2[j])
        s += '\n'

        for r in range(0, len(align.i1) + 2):
            j = 0
            s += '{} '.format(p1[r])

            for c in range(0, len(align.i2) + 2):
                s += "{} ".format(align.smat[r][c])
            s += "\n"

        s += '\npath\n\n'
        for r in range(0, len(align.i1) + 2):
            j = 0
            print()
            for c in range(0, len(align.i2) + 2):
                s += "{} ".format(align.pmat[r][c])
            s += "\n"

        return s

    def nw(self, ):
        """-----------------------------------------------------------------------------------------
        Real Needlman Wunsch algorithm with gaps as single steps. To get real NW without gap length
        penalty set self.gd = 0
        assumes that gap penalties gi and gd are negative

        s1 is horizontal, s2 is vertical. Origin is upper left (beginning of both sequences)
        :return:
        -----------------------------------------------------------------------------------------"""
        s1 = self.i1
        s2 = self.i2
        gd = self.gd
        gi = self.gi
        s = self.score.table
        not_possible = gi + max(len(s1), len(s2))

        # initialize
        previous = [0 for _ in range(len(s1))]
        current = [gi + j * gd for j in range(len(s1))]
        bestcol = [not_possible for _ in range(len(s1))]

        # previous row initialize with penalties to gap to [0,0]
        previous[0] = gi

        # main calculation
        for row in range(len(s2)):
            previous, current = current, previous

            # first position in row. first position in row can only make a gap in horizontal direction (bestrow)
            if row == 0:
                current[0] = s[s1[0]][s2[row]]
            else:
                current[0] = s[s1[0]][s2[row]] + gi

            bestrow = gi + gd * (row - 1)
            bestcol[0] = max( current[0] - gi, bestcol[0])

            for col in range(1, len(s1)):

                # find best score
                current[col] = s[s1[col]][s2[row]] + max(bestrow, bestcol[col - 1], previous[col - 1])

                # update bestrow for  next column in the current row
                if previous[col - 1] + gi > bestrow + gd:
                    bestrow = previous[col - 1] + gi
                else:
                    bestrow += gd

                # update bestcol for the next row
                if previous[col - 1] + gi > bestcol[col - 1] + gd:
                    bestcol[col - 1] = previous[col - 1] + gi
                else:
                    bestcol[col - 1] += gd

        return True

    def nwmod(self):
        """-----------------------------------------------------------------------------------------
        cur and paste from working version in main. not tested
        :return:
        -----------------------------------------------------------------------------------------"""
        # three pointer directions
        d = 1
        h = 2
        v = 4

        # score and path matrix are sequence length + 2
        align.smat = [[0 for c in range(len(align.i2) + 2)] for r in range(len(self.i1) + 2)]
        align.pmat = [[0 for c in range(len(align.i2) + 2)] for r in range(len(self.i1) + 2)]

        # fill first row (including end gap penalties)
        r = 0
        c = 0
        s = align.score.table
        i1 = align.i1
        i2 = align.i2
        gi = align.gi
        gd = align.gd
        smat = align.smat
        pmat = align.pmat

        # top and left edge conditions (end gaps)
        for c in range(1, len(align.i2) + 2):
            smat[r][c] = gi + c * gd
            pmat[r][c] = h
        c = 0
        for r in range(1, len(align.i1) + 2):
            smat[r][c] = gi + r * gd
            pmat[r][c] = v

        # body of comparison
        i = 0
        for r in range(1, len(align.i1) + 1):
            j = 0
            print()
            for c in range(1, len(align.i2) + 1):

                diag = smat[r - 1][c - 1] + s[i1[i]][i2[j]]
                left = smat[r][c - 1] + gd
                up = smat[r - 1][c] + gd
                best = max(diag, left, up)
                smat[r][c] = best
                print('r:{} c:{} i:{} j:{} s:{},{},{} diag:{} left:{} up:{} best:{}'.format(
                    r, c, i, j, i1[i], i2[j], s[i1[i]][i2[j]], diag, left, up, best))

                # set pointers
                if diag == best:
                    pmat[r][c] += d
                if left == best:
                    pmat[r][c] += h
                if up == best:
                    pmat[r][c] += v

                j += 1

            i += 1

        # final edge conditions (end gaps on right)
        r = len(align.i1) + 1
        end = len(align.i2)
        for c in range(1, len(align.i2) + 2):
            smat[r][c] = smat[r - 1][c - 1] + gi + (end - c + 1) * gd
            pmat[r][c] = d
        c = len(align.i2) + 1
        end = len(align.i1)
        for r in range(1, len(align.i1) + 2):
            smat[r][c] = smat[r - 1][c - 1] + gi + (end - r + 1) * gd
            pmat[r][c] = d

        print(align.matString())


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    align = Alignment()
    align.s1.seq = 'AGGC'
    # align.s1.seq = 'ACGTAAC'
    # align.s1.seq = "TAGATTTATCAT"
    print(align.s1.format())

    align.s2.seq = "AGCGT"
    # align.s2.seq = "CGAAGTC"
    # align.s2.seq = 'TATCATATGGT'
    print(align.s2.format())

    align.index()
    align.score.identity(pos=4, neg=-2)
    align.gi = -1
    align.gd = -1

    align.nw()

exit(0)

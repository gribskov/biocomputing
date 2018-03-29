class Alignment:
    """=============================================================================================
    simple dynamic programming without gaps
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        alignment constructor
        -----------------------------------------------------------------------------------------"""
        self.s1 = ''
        self.s2 = ''
        self.table = {}
        self.score = []
        self.pointer = []
        self.gapchar = '-'

    def globalNoGap(self):
        """-----------------------------------------------------------------------------------------
        Calculate score matrix for the no gap case
        :return: True
        -----------------------------------------------------------------------------------------"""
        # initialize score.  make matrix one position bigger than the sequence length.
        # since there are no gaps, the pointer matrix is never updated
        # s1 is horizontal, indexed by j (columns)
        # s2 is vertical, indexed by i (rows)
        self.score = [[0 for j in range(len(self.s1) + 1)] for i in range(len(self.s2) + 1)]
        self.pointer = [[0 for j in range(len(self.s1) + 1)] for i in range(len(self.s2) + 1)]

        # start at 1 to avoid special treatment for edges
        for i in range(1, len(self.score)):
            print('row j={}: {}'.format(i - 1, self.s2[1 - 1]))

            for j in range(1, len(self.score[i])):
                s = self.table[self.s2[i - 1]][self.s1[j - 1]]
                self.score[i][j] = self.score[i - 1][j - 1] + s
                print('ij:{},{} s={} S={}'.format(i, j, s, self.score[i][j]))

        return True

    def scoreMaxFinal(self):
        """-----------------------------------------------------------------------------------------
        retrun the maximum value and i, j position in the score matrix.  Only the final row and
        column position are checked by this method.  Note the score matrix has
        an extra row and column to deal with edge conditions
        :return: max, i, j
        -----------------------------------------------------------------------------------------"""
        imax = 1
        jmax = len(self.s1)
        maxscore = self.score[imax][jmax]

        # final column
        j = len(self.s1)
        for i in range(1, len(self.s2) + 1):
            if self.score[i][j] > maxscore:
                imax = i
                jmax = j
                maxscore = self.score[i][j]

        # final row
        i = len(self.s2)
        for j in range(1, len(self.s1) + 1):
            if self.score[i][j] > maxscore:
                imax = i
                jmax = j
                maxscore = self.score[i][j]

        return maxscore, imax - 1, jmax - 1

    def traceback(self, i, j):
        """-----------------------------------------------------------------------------------------
        Traceback the alignment from the specified position.  traceback is based on pointer matrix.
        :return: s1, s2, strings padded to align
        -----------------------------------------------------------------------------------------"""
        # construct strings backwards from starting position.  Overhanging sequences is filled with
        # gap characters.  i is the row index (s2), j is the column index (s1)

        # fill ends
        gap = self.gapchar
        als2 = ''
        als1 = ''
        ii = len(self.s2)
        while ii > i + 1:
            als2 += self.s2[ii - 1]
            als1 += gap
            ii -= 1

        jj = len(self.s1)
        while jj > j + 1:
            als1 += self.s1[jj - 1]
            als2 += gap
            jj -= 1

        while ii > 0 and jj > 0:
            als1 += self.s1[jj - 1]
            als2 += self.s2[ii - 1]

            if self.pointer[ii][jj] == 0:
                # no gap
                ii -= 1
                jj -= 1
            elif self.pointer[ii][jj] < 0:
                # vertical gap
                ii += self.pointer[ii][jj]
                jj -= 1
            else:
                # horizontal gap
                jj -= self.pointer[ii][jj]
                ii -= 1

        while ii > 0:
            als2 += self.s2[ii - 1]
            als1 += gap
            ii -= 1

        while jj > 0:
            als1 += self.s1[jj - 1]
            als2 += gap
            jj -= 1

        # reverse the strings when returning
        return als1[::-1], als2[::-1]


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    a = Alignment()
    a.s1 = 'AACGT'
    a.s2 = 'AACAGT'

    # scoring table with match=1 mismatch=-1
    a.table = {'A': {'A': 1, 'C': -1, 'G': -1, 'T': -1},
               'C': {'A': -1, 'C': 1, 'G': -1, 'T': -1},
               'G': {'A': -1, 'C': -1, 'G': 1, 'T': -1},
               'T': {'A': -1, 'C': -1, 'G': -1, 'T': 1}}

    a.globalNoGap()
    maxscore, i, j = a.scoreMaxFinal()
    print('maximum score: {}\t at {}, {}'.format(maxscore, i, j))

    al1, al2 = a.traceback(i, j)
    print(al1)
    print(al2)

    al1, al2 = a.traceback(5, 5)
    print(al1)
    print(al2)

    exit(0)

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

    def globalNoGap(self):
        """-----------------------------------------------------------------------------------------

        :return:
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

        return

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

        return maxscore, imax-1, jmax-1


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

    exit(0)

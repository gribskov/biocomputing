class Score:
    """=============================================================================================
    sequence scoring table object
    ============================================================================================="""
    import sys

    def __init__(self, alphabet='ACGT'):
        """-----------------------------------------------------------------------------------------
        score constructor
        -----------------------------------------------------------------------------------------"""
        self.alphabet = alphabet
        self.a2i = {}
        self.i2a = []
        self.table = []

        self.alphabetIndex()

    def alphabetIndex(self):
        """-----------------------------------------------------------------------------------------
        Initialze the a2i and i2f for converting between sequence characters and indices
        :return: alphabet size
        -----------------------------------------------------------------------------------------"""
        i = 0
        for char in self.alphabet:
            if char in self.a2i:
                # duplicate character
                sys.stderr.write(
                    'Score:alphabetIndex - duplicate alphabet character in {}\n'.format(alphabet))

            self.a2i[char] = i
            self.i2a.append(char)
            i += 1

        return i

    def identity(self, pos=1, neg=0):
        """-----------------------------------------------------------------------------------------
        define an identity matrix based on the current scoring matrix
        :return: True
        -----------------------------------------------------------------------------------------"""
        n = len(self.alphabet)
        self.table = [[neg for i in range(n)] for j in range(n)]
        for i in range(len(self.alphabet)):
            self.table[i][i] = pos

        return True

# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':

    score = Score()

    print('Default alphabet and indices')
    print('    alphabet:', score.alphabet)
    print('    a2i:', score.a2i)
    print('    i2a:', score.i2a)

    print('identity matrix')
    score.identity()
    print(score.table)

    exit(0)

import sys


class Score:
    """=============================================================================================
    sequence scoring table object
    ============================================================================================="""
    import sys

    def __init__(self, alphabet='ACGT'):
        """-----------------------------------------------------------------------------------------
        score constructor
        -----------------------------------------------------------------------------------------"""
        self.alphabetSet(alphabet=alphabet)

        # self.alphabet = alphabet
        # self.a2i = {}
        # self.i2a = []
        # self.table = []
        #
        # self.alphabetIndex()

    def __str__(self):
        return self.format(decimal=0)

    def format(self, space=2, decimal=0):
        """-----------------------------------------------------------------------------------------
        return a string with a formatted version of the scoring matrix
        :return: True
        -----------------------------------------------------------------------------------------"""
        str = ''

        # first find the width of the biggest value
        width = 0
        for i in range(len(self.alphabet)):
            for j in range(i, len(self.alphabet)):
                width = max(width, len(super.__str__(self.table[i][j])))

        width += space
        if decimal > 0:
            width += decimal + 1

        print('{:{width}}'.format('', width=width), end='')
        for i in self.alphabet:
            print('{0:>{width}}'.format(i, width=width), end='')
        print()

        for i in range(len(self.alphabet)):
            print('{:>{width}}'.format(self.i2a[i], width=width), end='')
            for j in range(len(self.alphabet)):
                print(
                    '{0:{width}.{decimal}f}'.format(self.table[i][j], width=width, decimal=decimal),
                    end='')
            print()

        return ''

    def alphabetIndex(self):
        """-----------------------------------------------------------------------------------------
        Initialze the a2i and i2f for converting between sequence characters and indices
        :return: alphabet size
        -----------------------------------------------------------------------------------------"""
        i = 0
        self.a2i = {}
        self.i2a = []
        for char in self.alphabet:
            if char in self.a2i:
                # duplicate character
                sys.stderr.write(
                    'Score:alphabetIndex - duplicate alphabet character in {}\n'.format(char))

            self.a2i[char] = i
            self.i2a.append(char)
            i += 1

        return i

    def alphabetSet(self, alphabet=''):
        """-----------------------------------------------------------------------------------------
        change the alphabet.
        1. rebulid a2i and i2f
        2. reset table to identity
        :return: alphabet size
        -----------------------------------------------------------------------------------------"""
        self.alphabet = alphabet
        self.alphabetIndex()
        self.identity()

        return len(self.alphabet)

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

    def strToIndex(self, str):
        """-----------------------------------------------------------------------------------------
        Convert a string to a list of indices in the current alphabet.  This speeds up score lookup

        :param str:
        :return: index list
        -----------------------------------------------------------------------------------------"""
        istring = []
        for c in str:
            istring.append(self.a2i[c])

        return istring

    def randomFromFrequency(self, freq):
        """-----------------------------------------------------------------------------------------
        Given a frequency probability vector, represented as a dict or list, construct a transition
        matrix that represents the probability of random transitions between letters.  Assuming the
        frequency vector is a row vector, this is the matrix with the transpose of the frequency
        vector as columns.

        :param freq: dict/list of float, frequencies of alphabet letters
        :return: True/False
        -----------------------------------------------------------------------------------------"""
        if isinstance(freq, dict):
            n = len(self.alphabet)
            for a in self.alphabet:
                i = self.a2i[a]
                for b in self.alphabet:
                    j = self.a2i[b]
                    self.table[i][j] = freq[a]
        elif isinstance(freq, list):
            for i in range(len(freq)):
                for j in range(len(freq)):
                    self.table[i][j] = freq[i]
        else:
            sys.stderr.write('Score:randomFromFrequency - unknown type of frequency vector\n')
            sys.stderr.write('frequency vector must be list or dict, found {}\n'.format(type(freq)))
            return False

        return True

    def conservedFromFrequency(self, freq, identity=1.0):
        """-----------------------------------------------------------------------------------------
        For a frequency probability vector and a conservation fraction construct a transition matrix
        in which identity fractionof letters are unchanged and 1-identity are randomly changed.

        This matrix is nether symmetric nor stationary

        :param freq: dict/list of float, frequencies of alphabet letters
        :param identity: float
        :return: True/False
        -----------------------------------------------------------------------------------------"""
        rfrac = (1.0 - identity)

        if isinstance(freq, dict):
            n = len(self.alphabet)

            for a in self.alphabet:
                i = self.a2i[a]
                self.table[i][i] = identity
                changed = 1.0 - freq[i]
                for b in self.alphabet:
                    j = self.a2i[b]
                    if i == j:
                        continue

                    self.table[j][i] = rfrac * freq[j] / changed

        elif isinstance(freq, list):
            for i in range(len(freq)):
                self.table[i][i] = identity
                changed = 1.0 - freq[i]
                for j in range(len(freq)):
                    if i == j:
                        continue

                    self.table[j][i] = rfrac * freq[j] / changed

        else:
            sys.stderr.write('Score:conservedFromFrequency - unknown type of frequency vector\n')
            sys.stderr.write('frequency vector must be list or dict, found {}\n'.format(type(freq)))
            return False

        return True

    def transition(self, freq):
        """-----------------------------------------------------------------------------------------
        Multiply the current transition matrix times the row vector, freq

        :param freq:
        :return: list, float; new frequencies
        -----------------------------------------------------------------------------------------"""
        new = [0 for _ in range(len(freq))]
        for i in range(len(freq)):
            for j in range(len(freq)):
                new[i] += freq[j] * self.table[j][i]

        return new


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

    print('print with one decimal')
    print(score.format(decimal=1))

    print('print in default str format')
    print(str(score))

    print('protein alphabet')
    asize = score.alphabetSet('ABCDEFGHIKLMNPQRSTVWXYZ')
    print('    alphabet size:', asize)
    print(str(score))

    print('convert sequence to index')
    sequence = 'ACEDFG'
    print('    sequence:', sequence)
    iseq = score.strToIndex(sequence)
    print('    iseq:', iseq)

    # transition matrices

    print('Transition matrix from dictionary')
    score = Score()
    freq = {'A': 0.4, 'C': 0.3, 'G': 0.2, 'T': 0.1}
    score.randomFromFrequency(freq)
    print(score.format(decimal=2))

    print('Transition matrix from list')
    score = Score()
    freq = [0.4, 0.3, 0.2, 0.1]
    score.randomFromFrequency(freq)
    print(score.format(decimal=2))

    # print('Transition from string - should fail')
    # score = Score()
    # freq = 'GATCGATC'
    # score.randomFromFrequency(freq)
    # print(score.format(decimal=2))

    print('90% conservation matrix')
    score = Score()
    freq = [0.4, 0.3, 0.2, 0.1]
    score.conservedFromFrequency(freq, 0.9)
    print(score.format(decimal=2))

    new = score.transition(freq)
    print('new frequencies: {}'.format(new))

    exit(0)

"""=================================================================================================
word composition of a string

Michael Gribskov     02 June 2021
================================================================================================="""
from itertools import product


class Word:

    def __init__(self, size=6, alphabet='ACGT', sequence=''):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.size = size
        self.initialized = False
        self.alphabet = alphabet
        self.sequence = sequence
        self.count = {}

    def counts_get(self):
        """-----------------------------------------------------------------------------------------
        count the words of length self.size, assumes word dict is initialized with all possible
        words, automatically initializes with zero if not previously done

        :return:int, number of words counted
        -----------------------------------------------------------------------------------------"""
        if not self.initialized:
            self.init_words()

        n = 0
        for pos in range(0, len(self.sequence) - self.size + 1, self.size):
            # print(self.sequence[pos: pos + self.size])
            word = self.sequence[pos: pos + self.size]
            self.count[word] += 1
            n += 1

        return n

    def counts_get_sparse(self):
        """-----------------------------------------------------------------------------------------
        count the words of length self.size. works with an uninitialized list of words
        :return:
        -----------------------------------------------------------------------------------------"""

        n = 0
        for pos in range(0, len(self.sequence) - self.size + 1, self.size):
            # print(self.sequence[pos: pos + self.size])
            word = self.sequence[pos: pos + self.size]
            if word in self.count:
                self.count[word] += 1
            else:
                self.count[word] = 1

            n += 1

        return n

    def init_words(self, initval=0):
        """-----------------------------------------------------------------------------------------
        initialize the word list as the number of permutations of the alphabet

        :param initval: int, initial value for counts, e.g., a1 for a plus one prior
        :return: int, total count of words
        -----------------------------------------------------------------------------------------"""
        n = len(self.alphabet) ** self.size * initval
        for w in product(self.alphabet, repeat=self.size):
            self.count[''.join(w)] = initval

        self.initialized = True

        return n


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    seq = 'GATGTCAAGGAAGATGAAGATGGACATGC'
    hexamer = Word(sequence=seq)
    c = hexamer.init_words(1)
    c += hexamer.counts_get()
    print(f'total count:{c}')
    print(hexamer.count)

    exit(0)

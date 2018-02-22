"""=================================================================================================
get a count of kmers in a sequence.  The real goal is to identify the infrequent kmers and use them
to probabilistically generate longer kmers unlikely to occur in the sequence.
================================================================================================="""
import sys
import itertools
from sequence.fasta import Fasta


class Kmer:
    """=============================================================================================
    Holds the kmer distribution

    ============================================================================================="""

    def __init__(self, k=8):
        """-----------------------------------------------------------------------------------------
        kmer constructor
        :param k: length of words
        -----------------------------------------------------------------------------------------"""
        self.k = k
        self.alphabet = 'ACGT'
        self.kmer = {}

    def setupWords(self, prior=1):
        """-----------------------------------------------------------------------------------------
        initialize all of the keys and give a count equal to prior
        :param prior:
        :return: nwords
        -----------------------------------------------------------------------------------------"""
        n = 0
        for word in itertools.product(self.alphabet, repeat=self.k):
            self.kmer[''.join(word)] = prior
            n += 1

        return n


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # fasta = Fasta()
    # fasta.open(sys.argv[1])

    KMER = 3
    kmer = Kmer(KMER)
    nwords = kmer.setupWords()
    print('{} words initialized'.format(nwords))

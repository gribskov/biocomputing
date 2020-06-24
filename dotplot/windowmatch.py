"""=================================================================================================


Michael Gribskov     24 June 2020
================================================================================================="""
from sequence.score import Score
from wordmatch import Base


class Windowmatch(Base, Score):
    """=============================================================================================
    Base provides some calculations related to diagonals and Fasta sequences
    Score provides scoring table methods

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Constructor for class windowmatch.  Sequences s1 and s2 are defined in Base class.

        -----------------------------------------------------------------------------------------"""
        Base.__init__(self)
        Score.__init__(self)
        self.window = 0
        self.threshold = 0

    def seqToInt(self):
        """-----------------------------------------------------------------------------------------
        Convert sequence strings to an integer arrays and stores in object.  An integer array is
        more convenient for direct lookups in the scoring table than a string

        :return: int, int length of sequence listss
        -----------------------------------------------------------------------------------------"""
        a2i = self.a2i

        self.i1 = [a2i[c] for c in self.s1.seq]
        self.i2 = [a2i[c] for c in self.s2.seq]

        return len(self.i1), len(self.i2)


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    from sequence.fasta import Fasta

    match = Windowmatch()
    print('done {}'.format(type(match)))
    print(match.alphabet)

    match.readNCBI('table/NUC4.4.matrix')
    print(match.format())

    fasta1 = Fasta(filename=sys.argv[1])
    fasta1.read()

    fasta2 = Fasta()
    fasta2.id = 'seq2'
    fasta2.doc = ' bases 1:50'
    fasta2.seq = fasta1.seq[:50]

    fasta1.seq = fasta1.seq[:200]

    match.s1 = fasta1
    match.s2 = fasta2
    l1, l2 = match.seqToInt()
    print(l1,l2)

    exit(0)

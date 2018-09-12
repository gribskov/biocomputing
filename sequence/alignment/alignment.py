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


    def index(self, i1=True, i2=True):
        """-----------------------------------------------------------------------------------------
        cretae an integer array version of the sequences.  this makes it easier to lookup scores.
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


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    align = Alignment()
    align.s1.seq = 'AGGC'
    print(align.s1.format())

    align.s2.seq = "AGCGT"
    print(align.s2.format())

    align.index()

    exit(0)

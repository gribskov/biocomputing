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


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    align = Alignment()
    align.s1.seq = 'aaa'
    print(align.s1.format())

    exit(0)

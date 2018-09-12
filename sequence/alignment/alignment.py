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
        self.gi = 0  # length independent gap penalty
        self.gd = 0  # length dependent gap penalty

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
    align.score.identity()
    align.gi = 0
    align.gd = -1

    smat = [[0 for c in range(len(align.i2))] for r in range(len(align.i1))]
    pmat = [[0 for c in range(len(align.i2))] for r in range(len(align.i1))]

    # fill first row (including end gap penalties)
    r = 0
    c = 0
    s = align.score.table
    i1 = align.i1
    i2 = align.i2
    gi = align.gi
    gd = align.gd

    smat[r][c] = s[i1[r]][i2[c]]
    for c in range(1, len(align.i2)):
        smat[r][c] = s[i1[r]][i2[c]] + gi + c * gd

    for r in range(1,len(align.i1)):
        c = 0
        smat[r][c] = s[i1[r]][i2[c]] + gi + r * gd
        for c in range(1, len(align.i2)):
            smat[r][c] = s[i1[r]][i2[c]]


exit(0)

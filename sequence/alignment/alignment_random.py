"""=================================================================================================
Significance of alignment using random simulation

Michael Gribskov     02 August 2020
================================================================================================="""
import sys
from sequence.fasta import Fasta
from sequence.score import Score

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fasta1 = Fasta(filename=sys.argv[1])
    fasta1.read()

    fasta2 = Fasta(filename=sys.argv[2])
    fasta2.read()

    cmp = Score()
    cmp.readNCBI('../../dotplot/table/BLOSUM62.matrix')

    exit(0)
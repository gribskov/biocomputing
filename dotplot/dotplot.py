"""=================================================================================================
dotplot

Michael Gribskov     17 June 2020
================================================================================================="""
import sys
from sequence.fasta import Fasta
from wordmatch import Match
from plotter import Plotter

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    
    # sequencename = sys.argv[1]
    # sys.stderr.write('\tID list: {}\n'.format(sequencename))
    # try:
    #     sequence = open(sequencename, 'r')
    # except:
    #     sys.stderr.write('\nUnable to open sequence ({})\n'.format(sequencename))
    #     exit(1)

    # read fasta formatted sequence
    fasta1 = Fasta(filename=sys.argv[1])
    fasta1.read()

    fasta2 = Fasta()
    fasta2.id = 'seq2'
    fasta2.doc = ' bases 1:50'
    fasta2.seq = fasta1.seq[:50]

    fasta1.seq = fasta1.seq[:200]

    match = Match()
    match.s1 = fasta1
    match.s2 = fasta2
    nmatch = match.identity()
    # fmatch = match.filterByLen(3)
    wmatch = match.filterByWindowCount(5,4)

    plot = Plotter()
    plot.match = match
    plot.setup()
    plot.lines(dither=0.03)
    # plot.dots()
    plot.show()


    exit(0)
"""=================================================================================================
Exam transcripts clustered at bundle (DN), component (c), gene (g), and isoform level to decide what
level best corresponds to a "gene"

Search is assumed to be diamond (blast tabular), tab separated
TRINITY_DN88428_c0_g1_i1   1146   1106   705   A0A059DJS1_EUCGR   290   1   135   135   65.2   4.2e-45   A0A059DJS1_EUCGR

qname qlen qbegin qend
sname slen sbegin send
alignlen score evalue stitle
================================================================================================="""
import sys
from blast import Blast

# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    infile = sys.argv[1]
    sys.stderr.write('Blast search: {}\n'.format(infile))
    blast = Blast(file=sys.argv[1])

    # fmt = 'qname qlen qbegin qend sname slen sbegin send alignlen score evalue stitle'
    fmt = 'qname sname id alignlen mismatch gapopen qbeg qend sbeg send evalue bit_score'
    nfields = blast.setFormat(fmt)

    exit(0)

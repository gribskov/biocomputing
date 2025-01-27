"""=================================================================================================
Exam transcripts clustered at bundle (DN), component (c), gene (g), and isoform level to decide what
level best corresponds to a "gene"

Search is assumed to be diamond (blast tabular), tab separated
See data/trinity_uniref_short.dmndblast
fields are:
qname qlen qbegin qend
sname slen sbegin send
alignlen score evalue stitle
================================================================================================="""
import sys
from blast import Blast

def readblock(blast, level, scores_query, skip=''):
    """---------------------------------------------------------------------------------------------
    read a group of sequences from the blast file that are the same at a certain level, b=bundle,
    c=component, g=gene, i=isoform

    :param blast: Blast object      Input search result with the format given in the header
    :param level: string            b, c, g, i
    :param skip: string             comma delimited list of keywords to skip
    :return: list                   parsed fields from matching lines
    ---------------------------------------------------------------------------------------------"""

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

    n = 0
    while blast.next():
        n += 1
        print('   ', n, blast.line)

    exit(0)

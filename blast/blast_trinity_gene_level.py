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
import re
from blast import Blast


class Group:
    """=============================================================================================
    group holds information about a set of trinity isoforms. the main items are a list of matching
    sequences (match) and keywords (keyword)

    Michael Gribskov
    ============================================================================================="""
    # re used for more complicated replacements (more than just a keyword)
    stitlere = re.compile(r'(n=\d+)')

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        group constructor
        -----------------------------------------------------------------------------------------"""
        self.id = ""
        self.level = ""
        self.match = []
        self.keyword = []
        self.stopword = ['UniRef90_', 'n=\d+', 'protein', 'uncharacterized', 'family',
                         'Tax=', 'TaxID=', 'RepID=']
        self.stopre = re.compile('|'.join(self.stopword))
        self.n = 0

    def isoform_add(self, blast, sid_prefix='UniRef90_'):
        """-----------------------------------------------------------------------------------------
        Adds the information from one line of the search result to the group

        :param blast: Blast object      contains one line of the search result
        :return: int                    number of items in the group
        -----------------------------------------------------------------------------------------"""
        self.id = splitID(blast.qid)
        self.level = 3
        self.match.append(blast.sid.replace(sid_prefix, ''))
        # filtered = self.stopre.sub('', blast.stitle)
        self.keyword.append(self.stopre.sub('', blast.stitle).split())
        self.n += 1

        return self.n

    @staticmethod
    def filter_keywords(string):
        """-----------------------------------------------------------------------------------------
        filter the stitle string
            remove first token (subject name)
            remove stopwords

        :param string: str      the hit information string from the search
        :return:
        -----------------------------------------------------------------------------------------"""
        pass


def readblock(blast, level, scores_query, skip=''):
    """---------------------------------------------------------------------------------------------
    read a group of sequences from the blast file that are the same at a certain level, b=bundle,
    c=component, g=gene, i=isoform

    :param blast: Blast object      Input search result with the format given in the header
    :param level: string            b, c, g, i
    :param skip: string             comma delimited list of keywords to skip
    :return: list                   parsed fields from matching lines
    ---------------------------------------------------------------------------------------------"""
    pass


# regex for splitID()
idre = re.compile(r'>*TRINITY_DN([^_]+)_c(\d+)_g(\d+)_i(\d+)')


def splitID(id):
    """---------------------------------------------------------------------------------------------
    Breakdown the trinity ID string to give the
    Cluster (bundle),  component, gene and isoform
    usage
        infohash = trinityID(trinity.id)

    :param id: string
    :return: dict       bundle, component, gene, isoform
    ---------------------------------------------------------------------------------------------"""
    cluster, component, gene, isoform = idre.match(id).groups()
    return {'bundle': cluster, 'component': component, 'gene': gene, 'isoform': isoform}


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    infile = sys.argv[1]
    sys.stderr.write('Blast search: {}\n'.format(infile))
    blast = Blast(file=sys.argv[1])

    fmt = 'qid qlen qbegin qend sid slen sbegin send alignlen pid score evalue stitle'
    # fmt = 'qname sname id alignlen mismatch gapopen qbeg qend sbeg send evalue bit_score'
    nfields = blast.setFormat(fmt)

    n = 0
    group = Group()
    while blast.next():
        n += 1
        print('   ', n, blast.line)
        group.isoform_add(blast)

    exit(0)

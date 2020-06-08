"""=================================================================================================
merge kraken2 reports into a table and filter by number of samples, counts of taxa present

 64.95  1895470 1895470 U       0       unclassified
 35.05  1022732 874     R       1       root
 34.34  1001967 29770   R1      131567    cellular organisms
 22.72  663149  34      D       2157        Archaea
 22.71  662861  95      P       28890         Euryarchaeota
 22.70  662418  236     P1      2290931         Stenosarchaea group
 22.67  661544  51822   C       183963            Halobacteria
 10.70  312220  13140   O       2235                Halobacteriales
  6.19  180549  9776    F       1963268               Haloarculaceae
  1.61  47094   1278    G       63743                   Natronomonas
  0.85  24773   0       S       416273                    Natronomonas moolapensis
  0.85  24773   24773   S1      268739                      Natronomonas moolapensis 8.8.11


columns:
1. % of minimizers mapping to clade
2. number of minimizers mapping to clade
3. number of minimizers mapping to taxon
4. taxonomic rank, note that integers may be appended to indicate subranks
    U unclassified
    R root
    D domain (kingdom)
    P phylum
    C class
    O order
    F family
    G genus
    S species
5. taxon ID
6. common name of taxon (indentation is significant)


Michael Gribskov     07 June 2020
================================================================================================="""
import sys


class Node:
    """=============================================================================================
    Single node in the multifurcating tree holding the taxonomy
    ============================================================================================="""

    # global variables to translate ranks to levels
    i2r = ['U', 'U1', 'R', 'R1', 'D', 'D1', 'P', 'P1', 'C', 'C1',
           'O', 'O1', 'F', 'F1', 'G', 'G1', 'S', 'S1']
    r2i = {}
    for i in range(len(i2r)):
        r2i[i2r[i]] = i

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.pct_mapped = 0.0
        self.n_mapped = 0
        self.n_taxon = 0
        self.rank = 'U'
        self.text = ''

        self.parent = []
        self.child = []

    @staticmethod
    def parse(line):
        """-----------------------------------------------------------------------------------------

        :param line: string, one line from taxonomy file
        :return: Node, new node in taxonomy
        -----------------------------------------------------------------------------------------"""
        (pct_mapped, n_mapped, n_taxon, rank, taxid, text) = line.rstrip().split('\t')

        node = Node()
        node.pct_mapped = pct_mapped
        node.n_mapped = n_mapped
        node.n_taxon = n_taxon
        node.rank = rank
        node.level = Node.r2i[rank]
        node.text = text

        return node


class Taxonomy():
    """=============================================================================================
    Single node in the multifurcating tree holding the taxonomy
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.index = {}
        self.root = None
        self.current = None

    def backtrack(self, node):
        """-----------------------------------------------------------------------------------------
        search backward through parents to find the first node with a smaller level (higher rank)
        This node will be the parent of the one at the provided level.
        Has no effect if level is greater than current.level

        :param level: integer, taxonomic level
        :return: Node, matching parent node
        -----------------------------------------------------------------------------------------"""
        while self.current.level >= node.level:
            self.current = self.current.parent

        return self.current


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    # open file
    # idlist
    reportname = sys.argv[1]
    sys.stderr.write('\tKraken2 report: {}\n'.format(reportname))
    try:
        report = open(reportname, 'r')
    except:
        sys.stderr.write('\nUnable to open Kraken2 report ({})\n'.format(reportname))
        exit(1)

    # define rank order for taxonomic ranks

    # read and store in taxonomy object
    tax = Taxonomy()

    # root
    line = report.readline()
    node = Node.parse(line)
    tax.root = node
    tax.current = node
    tax.index[node.text] = node

    for line in report:
        # create node and add to index
        node = Node.parse(line)
        tax.index[node.text] = node

        # set up parent child relationships. backtrack will set tax.current to the parent of the
        # new node

        tax.backtrack(node)
        node.parent = tax.current
        tax.current.child.append(node)
        tax.current = node

    exit(0)

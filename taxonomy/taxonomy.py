"""=================================================================================================
Two classes
    Node defines a single node in a taxonomy tree
    Taxonomy defines a tree and index

First implementation is for kraken2 format
Example, --report output of kraken2:

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
    K kingdom
    P phylum
    C class
    O order
    F family
    G genus
    S species
5. taxon ID
6. common name of taxon (indentation is significant)

Michael Gribskov     10 June 2020
================================================================================================="""
import sys


class Node:
    """=============================================================================================
    Single node in the multifurcating tree holding the taxonomy
    ============================================================================================="""

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
    def parseKraken(line):
        """-----------------------------------------------------------------------------------------
        parses one line of a Kraken taxonomy report; each line corresponds to a node in the taxonomy
        tree.  File should be tab delimited

            22.71  662861  95      P       28890         Euryarchaeota

        :param line: string, one line from taxonomy file
        :return: Node, new node in taxonomy
        -----------------------------------------------------------------------------------------"""
        (pct_mapped, n_mapped, n_taxon, rank, taxid, text) = line.rstrip().split('\t')

        node = Node()
        node.pct_mapped = pct_mapped
        node.n_mapped = n_mapped
        node.n_taxon = n_taxon
        node.rank = rank
        node.taxid = taxid
        node.text = text.lstrip()

        return node


class Taxonomy():
    """=============================================================================================
    The taxonomy tree
    ============================================================================================="""

    # global variables to translate ranks to numeric levels
    i2r = ['U', 'R', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S']
    r2i = {}
    for i in range(len(i2r)):
        r2i[i2r[i]] = i

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

        :param node: Node, child node for which a parent is sought
        :return: Node, matching parent node
        -----------------------------------------------------------------------------------------"""
        while self.rankGE(node):
            self.current = self.current.parent

        return self.current

    def dumpFromIndex(self):
        """-----------------------------------------------------------------------------------------
        Write out the tree in the order of the index.
        This should regenerate the original kraken2 output

        :return: integer, number of lines printed
        -----------------------------------------------------------------------------------------"""
        n = 0
        for taxon in self.index:
            node = self.index[taxon]
            # space = ' ' * node.level
            level = Taxonomy.r2i[node.rank[0]]
            space = ' ' * level
            print('{}\t{}\t{}\t{}\t{}{}'.format(node.pct_mapped, node.n_mapped,
                                                node.n_taxon, node.rank, space, node.text))
            n += 1

        return n

    def writeFormatted(self, min_percent, min_count, space=4, ):
        """-----------------------------------------------------------------------------------------
        Write out the tree in the order of the index.
        This should regenerate the original kraken2 output in better aligned format

        :param space: integer, number of spaces between columns
        :return: integer, number of lines
        printed
        -----------------------------------------------------------------------------------------"""
        # find the widths of the numeric columns
        width = {}
        width['n_mapped'] = 0
        width['n_taxon'] = 0
        for taxon in self.index:
            node = self.index[taxon]
            width['n_mapped'] = max(width['n_mapped'], len(node.n_mapped))
            width['n_taxon'] = max(width['n_taxon'], len(node.n_taxon))

        fmt = '{{:>{}}}\t{{:>{}}}\t{{:>{}}}\t{{:>{}}}\t{{}}{{}}'. \
            format(5 + space, width['n_mapped'] + space, width['n_taxon'] + space, 2 + space)

        n = 0
        for taxon in self.index:
            node = self.index[taxon]

            if min_count > int(node.n_taxon) and min_percent > float(node.pct_mapped):
                continue

            # space = '  ' * node.level
            level = Taxonomy.r2i[node.rank[0]]
            indent = ' ' * (space + level)

            print(fmt.format(node.pct_mapped, node.n_mapped,
                             node.n_taxon, node.rank, indent, node.text))
            n += 1

        return n

    # def rankLT(self, node):
    #     """-----------------------------------------------------------------------------------------
    #     Test if the current node has rank less than node
    #
    #     :param node: Node
    #     :return: Boolean
    #     -----------------------------------------------------------------------------------------"""
    #     crank = self.current.rank[0]
    #     cbaselevel = Taxonomy.r2i[crank]
    #
    #     nrank = node.rank[0]
    #     nbaselevel = Taxonomy.r2i[nrank]
    #
    #     if cbaselevel < nbaselevel:
    #         return True
    #
    #     elif cbaselevel == nbaselevel:
    #         if len(crank) < len(nrank):
    #             # only possibly when comparing a rank with suffix to one without
    #             return True
    #         elif len(nrank) == len(crank):
    #             # both are subranks of the same base level
    #             if int(crank[1:]) < int(nrank[1:]):
    #                 return True
    #
    #     return False

    def rankGE(self, node):
        """-----------------------------------------------------------------------------------------
        Test if the current node has rank greater or equal to node.

        :param node: Node
        :return: Boolean
        -----------------------------------------------------------------------------------------"""
        c0 = self.current.rank[0]
        c1 = self.current.rank[1:] or "0"
        # crank = self.current.rank[0]
        # clen = len(self.current.rank)
        # cbaselevel = Taxonomy.r2i[crank]

        n0 = node.rank[0]
        n1 = node.rank[1:] or "0"
        # nrank = node.rank[0]
        # nlen = len(node.rank)
        # nbaselevel = Taxonomy.r2i[nrank]

        if c0 == n0:
            if int(c1) < int(n1):
                return False
            else:
                return True

        else:
            return Taxonomy.r2i[n0] < Taxonomy.r2i[c0]

        # if self.current.rank == node.rank:
        #     return True
        #
        # if cbaselevel > nbaselevel:
        #     return True
        #
        # elif cbaselevel < nbaselevel:
        #     return False
        #
        # # both are subranks of the same base level
        # elif clen > nlen:
        #     # only possibly when comparing a rank with suffix to one without
        #     return True
        #
        # elif int(self.current.rank[1:]) >= int(node.rank[1:]):
        #
        #     return True

        return False

    def merge(self, tax):
        """-----------------------------------------------------------------------------------------
        Merge the taxonomy, tax, into this one.
        Method:
            Look up taxonomy names of tax in self.index
            If present:
                copy children from tax into self (no duplicates)
            else
                create a new node and add it as child of parent, and to index

        :param tax: Taxonomy
        :return: integer, number of nodes after merging
        -----------------------------------------------------------------------------------------"""

        # if this taxonomy is empty, copy the root node from tax
        if not self.index:
            self.root = Node()
            self.root.pct_mapped = tax.root.pct_mapped
            self.root.n_mapped = tax.root.n_mapped
            self.root.n_taxon = tax.root.n_taxon
            self.root.rank = tax.root.rank
            self.root.text = tax.root.text

            self.current = self.root
            self.index[self.root.text] = self.root

        for taxon in tax.index:
            node = tax.index[taxon]
            if taxon in self.index:
                # node exists in current taxonomy
                # print('{} exists'.format(taxon))
                new = self.index[taxon]
                new.pct_mapped += node.pct_mapped
                new.n_mapped += node.n_mapped
                new.n_taxon += node.n_taxon
                if not new.rank == node.rank:
                    sys.stderr.write('{}({}) : {}({})\t node ranks do not agree\n'. \
                                     format(new.text, new.rank, node.text, node.rank))
                # new.text = node.text
                if new.parent:
                    # avoids error when the root is reached
                    if not new.parent.taxid == node.parent.taxid:
                        sys.stderr.write('{}({}) : {}({})\t parents do not agree\n'. \
                                         format(new.text, new.parent.taxid,
                                                node.text, node.parent.taxid))
                    if new not in new.parent.child:
                        new.parent.child.append(new)

            else:
                # this node is not in the existing taxonomy so create it
                # sys.stderr.write('creating {} \n'.format(taxon))
                new = Node()
                new.pct_mapped = node.pct_mapped
                new.n_mapped = node.n_mapped
                new.n_taxon = node.n_taxon
                new.rank = node.rank
                new.text = node.text
                new.taxid = node.taxid
                if node.parent:
                    try:
                        new.parent = self.index[node.parent.taxid]
                    except:
                        sys.stderr.write('parent of {} ({}) is unknown\n'.format(
                            node.text, node.parent.taxid))
                    new.parent.child.append(new)
                self.index[node.taxid] = new
                # don't create children, this will be done when they are found

        return len(self.index)

    @staticmethod
    def parentRank(rank):
        """-----------------------------------------------------------------------------------------
        Apparently, there can be an unlimited number of sub levels to any taxonomic rank, indicated
        by numeric suffices: D > D1 > D2 ...
        This makes it nearly impossible to predfine a complete ordered list of levels. This method
        assumes if the rank is only one letter it can be looked up in r2i.  If the rank longer, the
        parent level will be the same rank letter with the suffix decremented by 1.  When the suffix
        is zero, it is dropped.

        :param rank: string, taxonomic rank
        :return: string, taxonmic rank of parent
        -----------------------------------------------------------------------------------------"""
        if len(rank) == 1:
            parent_rank = Taxonomy.i2r[Taxonomy.r2i[rank] - 1]
            return parent_rank

        # suffix = int(rank[1:]) - 1
        suffix = rank[1:]
        rank = rank[0]

        parent_rank = rank + suffix
        parent_rank.replace('0', '')

        return parent_rank

    @staticmethod
    def readKraken(report):
        """-----------------------------------------------------------------------------------------
        This constructor returns a taxonomy object based on a kraken report file

        See the header of this file  for format details

        :param file: filehandle, readable input file
        :return: Taxonomy object
        -----------------------------------------------------------------------------------------"""
        # first read the root and create the taxonomy
        tax = Taxonomy()

        line = report.readline()
        node = Node.parseKraken(line)
        tax.root = node
        tax.current = node
        tax.index[node.taxid] = node

        for line in report:
            # create node and add to index
            node = Node.parseKraken(line)
            tax.index[node.text] = node

            # set up parent child relationships. backtrack will set tax.current to the parent of the
            # new node

            tax.backtrack(node)
            node.parent = tax.current
            tax.current.child.append(node)
            tax.current = node

        return tax


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
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

    # read and store in taxonomy object
    tax = Taxonomy.readKraken(report)

    tax.writeFormatted(0.1, 1000, 2)
    exit(0)

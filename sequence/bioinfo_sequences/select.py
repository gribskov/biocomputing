"""=================================================================================================


Michael Gribskov     07 October 2025
================================================================================================="""
import sys


class Sequence():
    """=============================================================================================

    ============================================================================================="""

    def __init__(self, fh):
        """-----------------------------------------------------------------------------------------

        :param fh:
        -----------------------------------------------------------------------------------------"""
        self.fh = fh
        self.item = []

        self.readseqs()

    def readseqs(self):
        """-----------------------------------------------------------------------------------------
        read all sequences and store attributes. The attributes are stored as the string repr of a
        python dict so they can be converted with eval()

        :return:
        -----------------------------------------------------------------------------------------"""
        n_line = 0
        seq = None
        for line in self.fh:
            match n_line % 3:
                case 0:
                    # title line
                    seq = {'id': '', 'sequence': '', 'attr': []}
                    self.item.append(seq)
                    id, attrstr = line.split(' ', maxsplit=1)
                    seq['id'] = id.lstrip('>')
                    seq['attr'] = eval(attrstr)

                case 1:
                    # sequence
                    seq['sequence'] = line.rstrip()

                case 2:
                    # blank line
                    next

            n_line += 1

        return len(self.item)

    def index(self, key, nfield):
        """-----------------------------------------------------------------------------------------
        make an index use the attribute given as key. Use only the first nfield fields

        :param key: string      self.item[]['attr']['key']
        :param nfield: int      number of fields to use from the attribute, key
        :return:
        -----------------------------------------------------------------------------------------"""
        idx = {}
        n_item = 0
        for item in self.item:
            value = item['attr'][key]
            field = value.split()
            idxterm = field[nfield-1]
            if idxterm in idx:
                idx[idxterm].append(n_item)
            else:
                idx[idxterm] = [n_item]

            n_item += 1

        return idx


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    inpro = open(sys.argv[1], 'r')
    innuc = open(sys.argv[2], 'r')

    pro = Sequence(inpro)
    nuc = Sequence(innuc)

    genusn = nuc.index('organism_name', 1)
    genusp = pro.index('organism_name', 1)

    # sourcep = pro.index('pub_gene_id', 1)

    n_g = 0
    for g in genusp:
        n_g += 1

        pstr = '\t'
        if g in genusp:
            pstr = f'\t{genusp[g]}'
        nstr = '\t'
        if g in genusn:
            nstr = f'\t{genusn[g]}'

        print(f'{n_g}\t{g}\tp:{pstr}\tn:{nstr}')

    exit(0)

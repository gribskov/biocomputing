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

    def fasta(self, outfh, linelen=100):
        """-----------------------------------------------------------------------------------------
        write sequence to outfh in Fasta format

        :param outfh:
        :param linelen: int
        :return:
        -----------------------------------------------------------------------------------------"""


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    lenfrac = 0.1
    lower = 1.0 - lenfrac
    upper = 1.0 + lenfrac

    inpro = open(sys.argv[1], 'r')
    innuc = open(sys.argv[2], 'r')

    pro = Sequence(inpro)
    nuc = Sequence(innuc)

    genusn = nuc.index('organism_name', 1)
    genusp = pro.index('organism_name', 1)

    # average length
    n_pro = 0
    sum_len = 0
    for p in pro.item:
        sum_len += len(p['sequence'])
        n_pro += 1

    ave_len = sum_len / n_pro
    min_len = ave_len * lower
    max_len = ave_len * upper
    min_nuc = min_len * 3
    max_nuc = max_len * 3
    print(f'Average protein length: {ave_len:.1f}')


    # sourcep = pro.index('pub_gene_id', 1)

    n_g = 0
    for g in genusp:
        n_g += 1

        lenok = False
        for pid in genusp[g]:
            plen = len(pro.item[pid]['sequence'])
            if plen < min_len:
                # too short
                continue
            elif plen > max_len:
                # too long
                continue
            else:
                # length ok
                lenok = True

            if not lenok:
                # no suitable protein found, try next protein for this genus
                print(f'not suitable: {pid}\t{pro.item[pid]['attr']['organism_name']}\tlen:{len(pro.item[pid]['sequence'])}')
                continue

            # check the corresponding nucleotide sequence
            if pid < len(nuc.item):
                nlen = len(nuc.item[pid]['sequence'])
            else:
                # no matching nucleotide
                print(f'skipped/no nucleotide: {pid}\t{pro.item[pid]['attr']['organism_name']}\tlen:{len(pro.item[pid]['sequence'])}')
                break

            if nlen > min_nuc and nlen < max_nuc:
                # nucleotide length is ok
                # write out pro and nuc
                # skip other protein entries
                print(f'match {pid}\t{g}\t{pro.item[pid]['attr']['organism_name']}\tlen:',end=' ')
                print(f'{len(pro.item[pid]['sequence'])},{len(nuc.item[pid]['sequence'])}')
                break
            else:
                # nucleotide is not ok, try another protein
                continue

        # end of loop over proteins in genus



        # pstr = '\t'
        # if g in genusp:
        #     pstr = f'\t{genusp[g]}'
        # nstr = '\t'
        # if g in genusn:
        #     nstr = f'\t{genusn[g]}'
        #
        # print(f'{n_g}\t{g}\tp:{pstr}\tn:{nstr}')

    exit(0)

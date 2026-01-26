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

    def fasta(self, id, outfh, linelen=100, maxlen=1000000, long=True):
        """-----------------------------------------------------------------------------------------
        write sequence to outfh in Fasta format

        : param id: int         position of sequence in self.item
        :param outfh:           filehandle open for output
        :param linelen: int     length of lines for sequence
        :param maxlen: int      maximum sequence length (truncate at this length)
        :param long: bool       if true write the og group in the title
        :return: True
        -----------------------------------------------------------------------------------------"""
        seq = self.item[id]
        sequence = seq['sequence']
        name, rest = seq['attr']['organism_name'].split(' ', maxsplit=1)
        if long:
            outfh.write(f'>{name} {seq['attr']['og_name']} {seq['attr']['organism_name']}\n')
        else:
            outfh.write(f'>{name} {seq['attr']['organism_name']}\n')
        pos = 0
        seq_end = min(len(sequence), maxlen)
        while pos < seq_end:
            outfh.write(f'{sequence[pos:pos + linelen]}\n')
            pos += linelen

        return True


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    lenfrac = 0.05
    lower = 1.0 - lenfrac
    upper = 1.0 + lenfrac
    max_select = 200
    max_pro = 2000
    max_nuc = 6000

    inpro = open(sys.argv[1], 'r')
    innuc = open(sys.argv[2], 'r')

    pro = Sequence(inpro)
    nuc = Sequence(innuc)
    inpro.close()
    innuc.close()

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

    outpro = open('protein.out.fa', 'w')
    outnuc = open('nuc.out.fa', 'w')

    n_g = 0
    n_selected = 0
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
                n_selected += 1
                print(f'match {pid}\t{g}\t{pro.item[pid]['attr']['organism_name']}\tlen:',end=' ')
                print(f'{len(pro.item[pid]['sequence'])},{len(nuc.item[pid]['sequence'])}')
                pro.fasta(pid, outpro, maxlen=max_pro, long=False)
                nuc.fasta(pid, outnuc, maxlen=max_nuc, long=False)

                break

            else:
                # nucleotide is not ok, try another protein
                continue

            if n_selected >= max_select:
                break

        # end of loop over proteins in genus

    print(f'{n_selected} sequences selected')

    outpro.close()
    outnuc.close()

    exit(0)

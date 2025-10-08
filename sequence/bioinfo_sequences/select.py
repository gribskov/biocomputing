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




# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    inpro = open(sys.argv[1], 'r')
    innuc = open(sys.argv[2], 'r')

    pro = Sequence(inpro)
    nuc = Sequence(innuc)

    genus = nuc.index('organism_name')

    for p in pro.item:

    exit(0)

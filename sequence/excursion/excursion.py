"""-------------------------------------------------------------------------------------------------
Excursion analysis can identify regions that have a locally high density region of specific features

24 September 2018   Michael Gribskov

-------------------------------------------------------------------------------------------------"""
import sys
import re
from sequence.fasta import Fasta


class Feature:
    """---------------------------------------------------------------------------------------------
    Features to search for
    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        self.feature = '.'
        self.window = 0
        self.space = {}

    def count_space(self, fasta):
        """-----------------------------------------------------------------------------------------
        count the spacing between features.  There must be at least two occurences
        :param fasta:
        :return:
        -----------------------------------------------------------------------------------------"""
        result = re.finditer(r'(?=({}))'.format(self.feature), fasta)

        n = 0
        pos = None
        for hit in result:
            if pos:
                n += 1
                ln = hit.start(1) - pos
                print('ln', ln)
                try:
                    self.space['{}'.format(ln)] += 1
                except KeyError:
                    self.space['{}'.format(ln)] = 1

            else:
                pos = hit.start(1)

        return n


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    # fasta = Fasta(file=sys.argv[1])
    feature = Feature()
    feature.feature = "A.C"

    s = "AACGTACGTTGACAGCCCAACACGACGAGCATC"
    count = feature.count_space(s)
    print(count)

exit(0)

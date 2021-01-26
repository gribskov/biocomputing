"""=================================================================================================
Examin mpileup file produced from bam file by samtools

Michael Gribskov     26 January 2021
================================================================================================="""


class Mpileup:

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.file = ''
        self.fh = None
        self.parsed = {}

    # end of __init__

    def open_file(self, filename, mode):
        """-----------------------------------------------------------------------------------------
        Safe open for file

        :param filename: str, name of file
        :param mode:    file open mode, 'r', 'w', etc.
        :return: True if successful, False otherwise
        -----------------------------------------------------------------------------------------"""
        try:
            self.file = filename
            self.fh = open(filename, mode)
            return True

        except (OSError, IOError):
            print('mpileup::open - file open error ({}) in mode: {}'.format(filename, mode))
            return False

        # end of open_file

    def parse(self):
        """-----------------------------------------------------------------------------------------
        Read one line and return a hash of the tokens

        :return: dict, or False if read fails
        -----------------------------------------------------------------------------------------"""
        line = self.fh.readline()
        if line:
            token = line.split()
            # print(token)
            self.parsed = {'sequence': token[0],
                      'position': token[1],
                      'refbase': token[2],
                      'depth': token[3],
                      'bases': token[4],
                      'qualities': token[5]
                      }

            return self.parsed


        return False

    # end of parse


# End of Mpileup ===================================================================================

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    mp = Mpileup()
    mp.open_file('../data/PrFi_10k.mpileup', 'r')

    while mp.parse():
        print(mp.parsed['position'], mp.parsed['depth'])

    exit(0)

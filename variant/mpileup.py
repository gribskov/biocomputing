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
            self.fh = open(filename, mode )
            return True

        except (OSError, IOError):
            print('mpileup::open - file open error ({}) in mode: {}'.format(filename, mode))
            return False

        # end of open_file
# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    mp = Mpileup()
    mp.open_file('../data/PrFi_10k.mpileup', 'r')

    exit(0)

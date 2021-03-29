"""=================================================================================================
iterable fastq object

Michael Gribskov     29 March 2021
================================================================================================="""
import sys

class Fastq:

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------
        Object holds a single read plus the filename and filehandle (fh)
        -----------------------------------------------------------------------------------------"""
        self.id = ''
        self.sequence = ''
        self.quality = ''

        self.filename = filename
        self.fh = None

        if filename:
            try:
                self.fh = open(self.filename, 'r')
            except (IOError, OSError):
                sys.stderr.write(
                    'Fastq:__init__ - unable to open filename ({})\n'.format(self.filename))
                exit(1)

    def __iter__(self):
        return self

    def __next__(self):
        """-----------------------------------------------------------------------------------------
        Read the next 4 lines, a fastq read, from the filehandle, self.fh.

        :return: logical, True if an entry was read
        -----------------------------------------------------------------------------------------"""

        self.id = self.fh.readline().rstrip()
        self.sequence = self.fh.readline().rstrip()
        self.fh.readline().rstrip()  # skip separator line
        self.quality = self.fh.readline().rstrip()

        if self.quality:
            # if quality is None, eof was reached
            return True
        else:
            raise StopIteration
            return False


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    fastq = Fastq(filename='../../data/small.fq')

    n_read = 0
    for f in fastq:
        n_read += 1
        print('{}\t{}'.format(n_read, fastq.sequence))

    exit(0)

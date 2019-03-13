import sys


class Fastq:
    """=============================================================================================
    fastq.py

    A simple sequential fastq class
    ============================================================================================="""

    # class variable holds the GATCN alphabet
    complement = str.maketrans('acgtunACGTUN', 'tgcaanTGCAAN')

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------
        Object holds a single read plus the filename and filehandle (fh)
        -----------------------------------------------------------------------------------------"""
        self.id = ''
        self.sequence = ''
        self.separator = ''
        self.quality = ''
        self.error = '' # contains last error message

        self.filename = filename
        self.fh = None

        if filename:
            try:
                self.fh = open(self.filename, 'r')
            except (IOError, OSError):
                sys.stderr.write(
                    'Fastq:__init__ - unable to open filename ({})\n'.format(self.filename))
                exit(1)

    def next(self):
        """-----------------------------------------------------------------------------------------
        Read the next 4 lines, a fastq read, from the filehandle, self.fh.

        :return: logical, True if an entry was read
        -----------------------------------------------------------------------------------------"""

        self.id = self.fh.readline().rstrip()
        self.sequence = self.fh.readline().rstrip()
        self.separator = self.fh.readline().rstrip()
        self.quality = self.fh.readline().rstrip()

        if self.quality:
            # if quality is None, eof was reached
            return True
        else:
            return False

    def write(self, outfh):
        """-----------------------------------------------------------------------------------------
        Write out the current Fastq entry.  Formatting should be correct before calling this method

        :return: logical, True
       ------------------------------------------------------------------------------------------"""
        outfh.write('{}\n'.format(self.id))
        outfh.write('{}\n'.format(self.sequence))
        outfh.write('{}\n'.format(self.separator))
        outfh.write('{}\n'.format(self.quality))

        return True

    def sequence_revcomp(self):
        """-----------------------------------------------------------------------------------------
        reverse complement the sequence

        :return: string, reverse complement of sequence
        -----------------------------------------------------------------------------------------"""
        return self.sequence.translate(Fastq.complement)[::-1]

    def quality_rev(self):
        """-----------------------------------------------------------------------------------------
        Reverse the quality vector to correspond to reverse complement of sequence

        :return: list (reversed)
        -----------------------------------------------------------------------------------------"""
        return self.quality[::-1]

    def check(self):
        """-----------------------------------------------------------------------------------------
        Check the read format.  There are limited checks that can be done with fastq:
            ID line begins with @
            Sequence line has only ACGTN
            Separator line begins with +
            quality line is same length as sequence line (could check range of values)

        :return: logical, True if all checks are passed
        -----------------------------------------------------------------------------------------"""
        status = True
        self.error = ''

        if not self.id.startswith('@'):
            status = False
            self.error += 'ID line does not start with @\n'

        if not self.separator.startswith('+'):
            status = False
            self.error += 'Separator line does not start with +\n'

        if len(self.sequence) != len(self.quality):
            status = False
            self.error += 'Sequence and quality are different lengths +\n'

        tmp = self.sequence
        if len(tmp.strip('ACGTNacgtn')) != 0:
            status = False
            self.error += 'Sequence contains non ACGTN characters'

        return status


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    f1 = '../data/small.fq'
    fastq = Fastq(filename=f1)

    sys.stderr.write('\nBasic reading and writing:\n')
    n_read = 0
    while fastq.next():
        n_read += 1
        fastq.write(sys.stdout)

    sys.stderr.write('{} reads read from {}'.format(n_read, f1))

    sys.stderr.write('\nReverse first five reads:\n')
    n_read = 0
    fastq.fh.seek(0)
    while fastq.next():
        n_read += 1
        fastq.sequence = fastq.sequence_revcomp()
        fastq.quality = fastq.quality_rev()
        fastq.write(sys.stdout)
        if n_read >= 5:
            break

    sys.stderr.write('\ncheck for errors in reads:\n')
    while fastq.next():
        n_read += 1
        sys.stderr.write('\n\tno sequence\n')
        fastq.id = fastq.id.replace('@','X')
        if not fastq.check():
            sys.stdout.write(fastq.error)
            fastq.write(sys.stdout)

        if n_read >= 10:
            break

    exit(0)

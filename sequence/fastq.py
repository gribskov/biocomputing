class Fastq:
    """=============================================================================================
    fastq.py

    A simple sequential fastq class
    ============================================================================================="""
    import sys

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------
        Object holds a single read plus the filename and filehandle (fh)
        -----------------------------------------------------------------------------------------"""
        self.id = ''
        self.sequence = ''
        self.separator = ''
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

    def fastq_read(self):
        """-----------------------------------------------------------------------------------------
        Read the next 4 lines, a fastq read, from the filehandle, self.fh.

        :return: logical, True if an entry was read
        -----------------------------------------------------------------------------------------"""

        fq.id = fh.readline().rstrip()
        fq.sequence = fh.readline().rstrip()
        fq.separator = fh.readline().rstrip()
        fq.quality = fh.readline().rstrip()

        if fq.quality:
            # if quality is None, eof was reached
            return True
        else:
            return False

    def fastq_write(self):
        """-----------------------------------------------------------------------------------------
        Write out the current Fastq entry.  Formatting should be correct before calling this method

        :return: logical, True
       ------------------------------------------------------------------------------------------"""
        fh.write('{}\n'.format(self.id))
        fh.write('{}\n'.format(self.sequence))
        fh.write('{}\n'.format(self.separator))
        fh.write('{}\n'.format(self.quality))

        return True

    def sequence_revcomp(self):
        """-----------------------------------------------------------------------------------------
        reverse complement the sequence

        :return: string, reverse complement of sequence
        -----------------------------------------------------------------------------------------"""
        return self.sequence.translate(complement)[::-1]

    def quality_rev(selfl):
        """-----------------------------------------------------------------------------------------
        Reverse the quality vector to correspond to reverse complement of sequence

        :return: list (reversed)
        -----------------------------------------------------------------------------------------"""
        return self.quality[::-1]


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    exit(0)

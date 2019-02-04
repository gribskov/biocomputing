class Fastq():
    """=============================================================================================
    FastQ sequence data class. Data is assumed to be available in a file.

    Michael Gribskov     03 February 2019
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        FAstq object has attributes
        filename - filename of data file
        fh - filehandle of opend file
        id - ID of current entry
        sequence - sequence of current entry (string)
        quality - quality of current entry (list of int)
        -----------------------------------------------------------------------------------------"""
        self.filename = ''
        self.fh = None
        self.id = ''
        self.sequence = ''
        self.quality = []

    def openfile(self, fname):
        """-----------------------------------------------------------------------------------------
        Open the specified file fname.  Terminate with error status = 1 if there is an error.

        :param fname: string, filepath
        :return: filehandle
        -----------------------------------------------------------------------------------------"""
        self.filename = fname
        try:
            self.fh = open(fname, 'r')
        except (IOError, OSError):
            print('Error opening file {}'.format(fname))
            exit(1)

        return self.fh

    def next(self):
        """-----------------------------------------------------------------------------------------
        read the next entry from the file.  Return false whten all entries have been read

        :return: logical, false when done
        -----------------------------------------------------------------------------------------"""

        id = self.fh.readline()
        self.id = id.rstrip()

        seq = self.fh.readline()
        self.sequence = seq.rstrip()

        plus = self.fh.readline()
        qual = self.fh.readline()
        self.quality = self.get_quality_list(qual)

        return True

    def get_quality_list(self, qual):
        """---------------------------------------------------------------------------------------------
        Convert the quality string to a list of integer quality values.

        The ASCII characters in the quality line represent the probability that the base call at the
        corresponding position in the sequenceis incorrect. The are obtained as follows

        Q = int( -10 * log10( P(error) ) )

        Character = chr( 33 + Q )

        Q is referred to as the Phred quality, or more often as simply quality . The character '5',
        which is ASCII 53, represents quality = 20, which means an error probability of 10-2 or 0.01.
        A conventional choice for the minimum permissible quality is 20.

        Usage
            quality = get_quality_list(qual)

        :param qual: string
        :return: list of int
        ---------------------------------------------------------------------------------------------"""
        for qval in qual:
            q = ord(qval) - 33
            self.quality.append(q)

        return self.quality


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fastqname = '../HW1/8044.5k.fastq'
    fq = Fastq()

    fh = fq.openfile(fastqname)

    n_entry = 0
    while fq.next():
        n_entry += 1

        print('\t{}\t{}'.format(n_entry, fq.id))
        print('\t{}'.format(fq.sequence))

        if n_entry > 4:
            break

    # end of loop of Fastq entries

    exit(0)

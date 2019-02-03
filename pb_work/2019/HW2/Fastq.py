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


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fastqname = '../HW1/8044.5k.fastq'
    fq = Fastq()

    fh = fq.openfile(fastqname)
    print('file: {}'.format(fq.filename))
    print('filehandle: {}'.format(fh))

    exit(0)

"""=================================================================================================
Calculate mutations per read from the CIGAR string and MD: field in a sam file

Michael Gribskov     09 March 2021
================================================================================================="""
import sys


class Sam:
    """---------------------------------------------------------------------------------------------

    ---------------------------------------------------------------------------------------------"""

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------
        Define attributes for the 11 mandatory fields.  Anything else is stored in the optional
        attribute
        -----------------------------------------------------------------------------------------"""
        self.filename = ''
        self.fh = None

        self.qname = ''
        self.flag = 0
        self.rname = ''
        self.pos = None
        self.mapq = None
        self.cigar = ''
        self.rnext = ''
        self.pnext = None
        self.tlen = None
        self.seq = ''
        self.qual = ''
        self.option = {}

        if filename:
            self.filename = filename
            try:
                self.fh = open(filename, 'r')
            except (OSError, IOError):
                sys.stderr.write('Sam::__init__ - error opening SAM file ({})'.format(filename))

    def parse_alignment(self):
        """-----------------------------------------------------------------------------------------
        Parse the tab delimited SAM format.  Mandatory fields are stored inattributes, optional
        fields are stored in the self.option attribute which is a dictionary

        :return:
        -----------------------------------------------------------------------------------------"""
        line = self.fh.readline()
        if not line:
            return False

        field = line.rstrip().split('\t')
        self.qname = field[0]
        self.flag = int(field[1])
        self.rname = field[2]
        self.pos = int(field[3])
        self.mapq = int(field[4])
        self.cigar = field[5]
        self.rnext = field[6]
        self.pnext = int(field[7])
        self.tlen = int(field[8])
        self.seq = field[9]
        self.qual = field[10]
        self.option = {}

        for i in range(11, len(field)):
            tag, type, value = field[i].split(':')
            self.option[tag] = [type, value]

        return True


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    sam = Sam(filename='chr11.test.sam')
    line = sam.fh.readline()
    while line.startswith('@'):
        line = sam.fh.readline()

    n = 0
    while n < 10:
        sam.parse_alignment()

    n += 1

    exit(0)

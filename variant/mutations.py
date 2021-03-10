"""=================================================================================================
Calculate mutations per read from the CIGAR string and MD: field in a sam file

Michael Gribskov     09 March 2021
================================================================================================="""
import sys
import re


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

        self.header = {}

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

    def read_header(self):
        """-----------------------------------------------------------------------------------------
        Doesn't do much with the header, just parses out whatever tags are there and saves them
        in lists

        :return:
        -----------------------------------------------------------------------------------------"""
        line = sam.fh.readline().rstrip()
        while line.startswith('@'):
            field = line.split('\t')
            tag = field[0].replace('@','')
            if tag in self.header:
                self.header[tag].append(field[1:])
            else:
                self.header[tag] = [field[1:]]

            line = sam.fh.readline().rstrip()
            self.line = line

        return

    def parse_alignment(self):
        """-----------------------------------------------------------------------------------------
        Parse the tab delimited SAM format.  Mandatory fields are stored inattributes, optional
        fields are stored in the self.option attribute which is a dictionary

        :return:
        -----------------------------------------------------------------------------------------"""
        if not self.line:
            return False

        field = self.line.rstrip().split('\t')
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

        self.line = self.fh.readline()

        return True

    def count_diff(self):
        """-----------------------------------------------------------------------------------------
        count the number of difference between this read and reh reference.  The count is the number
        of alphabetic characters in the MD string, plus the number of insertions in the CIGAR
        string. Deletions in the CIGAR string appear in the MD string so they can be ignored.

        :return: int, number of sequence differences
        -----------------------------------------------------------------------------------------"""
        diff = 0

        insert = re.compile('(\d+)I')
        found = insert.findall(self.cigar)
        for f in found:
            diff += int(f)

        md = self.option['MD'][1]
        for base in 'ACGTN':
            diff += md.count(base)

        return diff

        return

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    sam = Sam(filename='chr11.test.sam')
    sam.read_header()

    n = 0
    while sam.parse_alignment():
        diff = sam.count_diff()
        print('{}\t{}\t{}\t{}\t{}'.format(sam.qname, sam.pos, sam.cigar, sam.option['MD'][1], diff))
        n += 1
        if n > 30:
            break

    exit(0)

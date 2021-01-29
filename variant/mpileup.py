"""=================================================================================================
Examin mpileup file produced from bam file by samtools

Michael Gribskov     26 January 2021
================================================================================================="""
import re
from math import log


class Mpileup:

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.file = ''
        self.fh = None
        self.parsed = {}
        self.count = {}

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
                           'position': int(token[1]),
                           'refbase': token[2],
                           'depth': int(token[3]),
                           'bases': token[4],
                           'qualities': token[5]
                           }

            return self.parsed

        return False

    # end of parse

    def countchar(self):
        """-----------------------------------------------------------------------------------------
        count the sequence characters at a position, forward and reverse strand, read ends, etc
        :return:
        -----------------------------------------------------------------------------------------"""
        readstart = re.compile(r'\^(\S)')
        indel = re.compile(r'([+-])([0-9]+)([ACGTNacgtn*#]+)')
        # indels = re.compile(r'[+-]([0-9]+)[ACGTNacgtn*#]{\g<1>}')

        bases = self.parsed['bases']
        count = {'indel': [], 'mapqual': [],
                 'A': 0, 'a': 0, 'C': 0, 'c': 0, 'G': 0, 'g': 0, 'T': 0, 't': 0, 'N': 0, 'n': 0,
                 '*': 0, '$': 0,
                 'forward': 0, 'backward': 0, 'strand_bias': 0}
        indel_list = []

        # first check for multi-character sequences
        begin = 0
        if '^' in bases:
            # beginning of sequence marker and map quality
            for match in readstart.finditer(bases):
                # print(bases)
                # print(match.group(1), ord(match.group(1)) - 33)
                bases = readstart.sub('', bases)
                begin += 1

        if '-' in bases or '+' in bases:
            # indel present
            match = indel.search(bases)
            while match:
                mlen = int(match.group(2))
                indel_list.append(mlen)
                old = '{}{}{}'.format(match.group(1), mlen, match.group(3)[:mlen])
                # print(bases, old)
                bases = bases.replace(old, '', 1)
                match = indel.search(bases)

            # print('{}\t{}'.format(self.parsed['depth'], bases))

        for c in bases:
            if c in count:
                count[c] += 1
            else:
                count[c] = 1

        # count the upper cas (forward strand) and lower case (backward strand)
        count['forward'] = 0
        count['backward'] = 0
        for c in 'acgtnACGTN':
            if c in count:
                if c.islower():
                    count['forward'] += count[c]
                elif c.isupper():
                    count['backward'] += count[c]

        # bias with +1 prior
        count['strand_bias'] = log((count['forward'] + 1) / (count['backward'] + 1), 2)
        count['end'] = begin + count['$']
        count['end_frac'] = (count['end'] + 1) / (self.parsed['depth'] + 1)

        count['indel'] = len(indel_list) + count['*']
        count['indel_frac'] = (count['indel'] + 1) / (self.parsed['depth'] + 1)
        self.count = count
        return self.count

    # end of count

    def mendel(self, hommax=0.1, hetmin=0.4, hetmax=0.6):
        """-----------------------------------------------------------------------------------------
        Decides if a position is homozygous (H) or heterozygous (h).  Positions with less than
        hommax alternate bases are H, positions with between hetmin and hetmax alternate bases are
        h, anything els is N (non-mendelian).  If no base is present, an empty string is returned.

        :param hommax: float, maximum alternate allele frequency for homozygous
        :param hetmin: float, minimum alternate allele frequency for heterozygous
        :param hetmax: float, maximum alternate allele frequency for heterozygous
        :return: string, H, h, N, or ''
        -----------------------------------------------------------------------------------------"""
        genotype = ''
        depth = self.parsed['depth']

        if depth == 0:
            # no base present
            return genotype

        # find the maximum allele
        depth += 1  # plus 1 prior
        maxfreq = 0
        maxallele = ''
        genotype = 'N'
        for c in 'ACGTN':
            freq = (count[c] + count[c.lower()] + 1) / depth
            if freq > maxfreq:
                maxfreq = freq
                maxallele = c
                minorfreq = 1 - freq
                if minorfreq < hommax:
                    genotype = 'H'
                elif minorfreq > hetmin and minorfreq < hetmax:
                    genotype = 'h'

        return genotype

    # end of mendel


# End of Mpileup ===================================================================================

# --------------------------------------------------------------------------------------------------
#

# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    depth_min = 10
    bias_max = 3

    mp = Mpileup()
    mp.open_file('../data/PrFi_10k.mpileup', 'r')

    total = {'indel': []}

    while mp.parse():
        count = mp.countchar()
        # print('\t', count)
        genotype = mp.mendel()

        status = ''
        if mp.parsed['depth'] < depth_min:
            status += 'L'
        if abs(count['strand_bias']) >= bias_max:
            status += 'B'

        print('{}\t{}\t{}:{}:{}:{}:{}\t{}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{}'.format(
            mp.parsed['position'],
            mp.parsed['depth'],
            count['A'] + count['a'],
            count['C'] + count['c'],
            count['G'] + count['g'],
            count['T'] + count['t'],
            count['N'] + count['n'],
            genotype,
            count['strand_bias'],
            count['end_frac'],
            count['indel'],
            count['indel_frac'],
            status
        ))

    exit(0)

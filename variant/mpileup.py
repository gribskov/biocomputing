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
                 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, '*': 0, '$': 0}

        # first check for multi-character sequences
        begin = 0
        if '^' in bases:
            # marks beginning of sequence and map shows quality
            for match in readstart.finditer(bases):
                # print(bases)
                # print(match.group(1), ord(match.group(1)) - 33)
                bases = readstart.sub('', bases)
                begin += 1

        indel_list = []
        if '-' in bases or '+' in bases:
            # indel present
            match = indel.search(bases)
            while match:
                mlen = int(match.group(2))
                if match.group(1) == '+':
                    # only add insertions, the deletions are marked as *
                    indel_list.append(mlen)
                old = '{}{}{}'.format(match.group(1), mlen, match.group(3)[:mlen])
                # print(bases, old)
                bases = bases.replace(old, '', 1)
                match = indel.search(bases)

            # print('{}\t{}'.format(self.parsed['depth'], bases))

        # count the bases and reflect into upper case. we only need the forward backward info for
        # the aggregate bases
        count['forward'] = 0
        count['backward'] = 0
        for c in bases:
            # each base at this position in the genome: should be acgtnACGTN*$
            if c.islower():
                count['backward'] += 1
                count[c.upper()] += 1
            else:
                count['forward'] += 1
                count[c] += 1

        # # count the upper case (forward strand) and lower case (backward strand)
        # for c in 'ACGTN':
        #     if c in count:
        #         if c.islower():
        #             count['forward'] += count[c]
        #         elif c.isupper():
        #             count['backward'] += count[c.lower()]

        # bias with +1 prior
        count['strand_bias'] = log((count['forward'] + 1) / (count['backward'] + 1), 2)
        count['end'] = begin + count['$']
        count['end_frac'] = (count['end'] + 1) / (self.parsed['depth'] + 1)

        count['indel'] = len(indel_list) + count['*']
        count['indel_frac'] = (count['indel'] + 1) / (self.parsed['depth'] + 1)
        self.count = count
        return self.count

    # end of count

    def mendel(self, hommax=0.12, hetmin=0.4, hetmax=0.6, known_min=0.89):
        """-----------------------------------------------------------------------------------------
        Decides if a position is homozygous (H) or heterozygous (h).  Positions with less than
        hommax alternate bases are H, positions with between hetmin and hetmax alternate bases are
        h, anything els is N (non-mendelian).  If no base is present, an empty string is returned.

        :param hommax: float, maximum alternate allele frequency for homozygous
        :param hetmin: float, minimum alternate allele frequency for heterozygous
        :param hetmax: float, maximum alternate allele frequency for heterozygous
        :param known_min: float, minimum fraction of known bases
        :return: string, H, h, N, or ''
        -----------------------------------------------------------------------------------------"""
        genotype = 'U'
        major = '.'
        minor = '.'
        depth = self.parsed['depth']
        freq = {}

        if depth == 0:
            # no base present
            return genotype, major, minor

        # find the maximum allele
        genotype = 'N'
        depth += 1  # plus 1 prior
        for c in 'ACGTN*':
            freq[c] = (count[c] + 0.167) / depth

        freq_order = sorted(freq.keys(), key=lambda x: freq[x], reverse=True)
        major = freq_order[0]
        minor = freq_order[1]

        if freq[major] + freq[minor] < known_min:
            genotype = 'U'
        elif 1 - freq[major] < hommax:
            genotype = 'H'
        elif freq[minor] > hetmin and freq[minor] < hetmax:
            genotype = 'h'
        else:
            genotype = 'N'

        return genotype, major, minor

    # end of mendel


# End of Mpileup ===================================================================================

# --------------------------------------------------------------------------------------------------
#

# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    depth_min = 10
    bias_max = 3
    indel_max = 0.03

    mp = Mpileup()
    mp.open_file('../data/PrFi_10k.mpileup', 'r')

    total = {'indel': []}

    while mp.parse():
        count = mp.countchar()
        # print('\t', count)
        genotype, major, minor = mp.mendel()
        if genotype and 'H' in genotype:
            continue

        status = ''
        if mp.parsed['depth'] < depth_min:
            status += 'L'
        if abs(count['strand_bias']) >= bias_max:
            status += 'B'
        if count['indel_frac'] > indel_max:
            status += 'I'

        print('{}\t{}\t{}:{}:{}:{}:{}:{}\t{} {}{}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{}'.format(
            mp.parsed['position'],
            mp.parsed['depth'],
            count['A'], count['C'], count['G'], count['T'], count['N'], count['*'],
            genotype, major, minor,
            count['strand_bias'],
            count['end_frac'],
            count['indel'], count['indel_frac'],
            status
        ))

    exit(0)

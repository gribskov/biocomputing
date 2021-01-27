"""=================================================================================================
Examin mpileup file produced from bam file by samtools

Michael Gribskov     26 January 2021
================================================================================================="""
import re


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
        count = {'indel': [], 'mapqual': []}
        indel_list = []

        # first check for multi-character sequences
        if '^' in bases:
            # beginning of sequence marker and map quality
            for match in readstart.finditer(bases):
                # print(bases)
                # print(match.group(1), ord(match.group(1)) - 33)
                bases = readstart.sub('', bases)

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

        # error checking
        #     if self.parsed['depth'] != len(bases):
        #     dollar = bases.count('$')
        #     if len(bases) - dollar != self.parsed['depth']:
        #         print('{}:{}\t{}\t{}'.format(self.parsed['position'],
        #                                      self.parsed['depth'],
        #                                      len(bases),
        #                                      bases))
        count['indel'] = indel_list
        self.count = count
        return self.count

    # end of count


# End of Mpileup ===================================================================================

# --------------------------------------------------------------------------------------------------
#

# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    mp = Mpileup()
    mp.open_file('../data/PrFi_10k.mpileup', 'r')

    total = {'indel': []}
    while mp.parse():
        print('=>', mp.parsed['position'], mp.parsed['depth'])
        # print(mp.countchar())
        mp.countchar()
        for key in mp.count:
            if key is 'indel':
                if mp.count['indel']:
                    total['indel'] += mp.count['indel']

            elif key in total:
                total[key] += mp.count[key]

            else:
                total[key] = mp.count[key]

    print(total)

    exit(0)

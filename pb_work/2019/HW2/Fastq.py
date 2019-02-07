class Fastq():
    """=============================================================================================
    FastQ sequence data class. Data is assumed to be available in a file.

    Michael Gribskov     03 February 2019
    ============================================================================================="""

    def __init__(self, filename = ''):
        """-----------------------------------------------------------------------------------------
        FAstq object has attributes
        filename - filename of data file
        fh - filehandle of opened file
        id - ID of current entry
        sequence - sequence of current entry (string)
        quality - quality of current entry (list of int)
        -----------------------------------------------------------------------------------------"""
        self.filename = filename
        self.fh = None
        self.id = ''
        self.sequence = ''
        self.quality = []
        self.quality_min = 0

        if filename:
            self.openfile()

    def openfile(self):
        """-----------------------------------------------------------------------------------------
        Open the specified file fname.  Terminate with error status = 1 if there is an error.

        :param fname: string, filepath
        :return: filehandle
        -----------------------------------------------------------------------------------------"""
        try:
            self.fh = open(self.filename, 'r')
        except (IOError, OSError):
            print('Error opening file ({})'.format(self.filename))
            exit(1)

        return self.fh

    def next(self):
        """-----------------------------------------------------------------------------------------
        read the next entry from the file.  Return false when all entries have been read

        :return: logical, false when done
        -----------------------------------------------------------------------------------------"""

        id = self.fh.readline()
        self.id = id.rstrip()

        seq = self.fh.readline()
        self.sequence = seq.rstrip()

        plus = self.fh.readline()
        qual = self.fh.readline()
        if not qual:
            # at end of file nothing can be read, qual is empty and therefore False
            return False

        self.quality = self.get_quality_list(qual)

        return True

    def get_quality_list(self, qual):
        """-----------------------------------------------------------------------------------------
        Convert the quality string to a list of integer quality values.

        The ASCII characters in the quality line represent the probability that the base call at the
        corresponding position in the sequence is incorrect. The are obtained as follows

        Q = int( -10 * log10( P(error) ) )

        Character = chr( 33 + Q )

        Q is referred to as the Phred quality, or more often as simply quality . The character '5',
        which is ASCII 53, represents quality = 20, which means an error probability of 10-2 or 0.01.
        A conventional choice for the minimum permissible quality is 20.

        Usage
            quality = get_quality_list(qual)

        :param qual: string
        :return: list of int
        -----------------------------------------------------------------------------------------"""
        self.quality = []
        for qval in qual:
            q = ord(qval) - 33
            self.quality.append(q)

        return self.quality

    def base_count(self):
        """-----------------------------------------------------------------------------------------
        Count the bases in the sequence that have quality > self.quality_min

        :return: dict of int, keys are bases plus 'total'
        -----------------------------------------------------------------------------------------"""
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'total': 0}
        for pos in range(len(self.sequence)):
            if self.quality[pos] >= self.quality_min:
                try:
                    count[self.sequence[pos]] += 1
                    count['total'] += 1
                except KeyError:
                    count[self.sequence[pos]] = 1
                    count['total'] = 1

        return count

    def trim(self):
        """-----------------------------------------------------------------------------------------
        Truncate the sequence and quality before the first base with quality < quality_min. This
        operation cannot be undone.

        :return: int, trimmed length
        -----------------------------------------------------------------------------------------"""
        pos = 0
        for q in self.quality:
            if q < self.quality_min:
                break
            pos += 1

        self.sequence = self.sequence[:pos]
        self.quality = self.quality[:pos]

        return pos


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fastqname = '../HW1/8044.5k.fastq'
    fq = Fastq(fastqname)

    n_entry = 0
    all = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'total': 0}
    hq = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'total': 0}
    trim = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'total': 0}
    while fq.next():
        n_entry += 1

        # print(n_entry, fq.sequence)
        fq.quality_min = 0
        count = fq.base_count()
        for base in count:
            all[base] += count[base]

        fq.quality_min = 20
        count = fq.base_count()
        for base in count:
            hq[base] += count[base]

        fq.trim()
        count = fq.base_count()
        for base in count:
            trim[base] += count[base]

        # if n_entry > 400:
        #     break

    # end of loop of Fastq entries
    print('{} entries read'.format(n_entry))
    print('\nAll bases')
    for base in all:
        print('{:>10s}: {}'.format(base, all[base]))

    print('\nHigh quality bases')
    for base in hq:
        print('{:>10s}: {}'.format(base, hq[base]))

    print('\nTrimmed bases')
    for base in hq:
        print('{:>10s}: {}'.format(base, trim[base]))

    exit(0)

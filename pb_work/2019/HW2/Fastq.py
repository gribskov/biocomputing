class Fastq():
    """=============================================================================================
    FastQ sequence data class. Data is assumed to be available in a file.

    Michael Gribskov     03 February 2019
    ============================================================================================="""

    def __init__(self, filename=''):
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
        Open the specified file in self.filename.  Terminate with error status = 1 if there is an
        error.

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


class Count:
    """=============================================================================================
    The count object keeps a count of letters, generally from a DNA or protein sequence

    ============================================================================================="""

    def __init__(self, alphabet='ACGTN'):
        """-----------------------------------------------------------------------------------------
        Constructor
        Letter counts are stored in a dictionary  with the letters as keys
        The alphabet specifies the valid letters that can be counted
        -----------------------------------------------------------------------------------------"""
        self.alphabet = alphabet
        self.count = {}
        for letter in alphabet:
            self.count[letter] = 0
        self.total = 0

    def add(self, count):
        """-----------------------------------------------------------------------------------------
        Add the counts in count to the object

        :param count: dict
        :return: int, total counts
        -----------------------------------------------------------------------------------------"""
        total = 0
        for letter in self.alphabet:
            self.count[letter] += count[letter]
            total += count[letter]

        self.total += total
        return total

    def report(self, title=''):
        """-----------------------------------------------------------------------------------------
        Return a formatted string reporting the counts

        :param title: string, title to print before counts
        :return: string
        -----------------------------------------------------------------------------------------"""
        report = '\n{}\n'.format(title)
        total = 0
        for base in self.alphabet:
            report += '{:>10s}: {}\n'.format(base, self.count[base])
            total += self.count[base]

        report += '{:>10s}: {}\n'.format('total', total)

        return report


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fastqname = '../HW1/8044.5k.fastq'
    fq = Fastq(fastqname)

    n_entry = 0
    all_bases = Count()   # all bases
    hq_bases = Count()    # bases with quality >= minimum
    trimmed_bases = Count()  # bases trimmed at firt base with quality < minimum
    while fq.next():
        n_entry += 1

        fq.quality_min = 0
        all_bases.add(fq.base_count())

        fq.quality_min = 20
        hq_bases.add(fq.base_count())

        fq.trim()
        trimmed_bases.add(fq.base_count())

    # end of loop of Fastq entries
    print('{} entries read'.format(n_entry))
    print(all_bases.report('All bases'))
    print(hq_bases.report('High quality bases'))
    print(trimmed_bases.report('Trimmed bases'))

    exit(0)

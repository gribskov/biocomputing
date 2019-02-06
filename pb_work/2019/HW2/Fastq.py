class Fastq():
    """=============================================================================================
    FastQ sequence data class. Data is assumed to be available in a file.

    Michael Gribskov     03 February 2019
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        FAstq object has attributes
        filename - filename of data file
        fh - filehandle of opened file
        id - ID of current entry
        sequence - sequence of current entry (string)
        quality - quality of current entry (list of int)
        -----------------------------------------------------------------------------------------"""
        self.filename = ''
        self.fh = None
        self.id = ''
        self.sequence = ''
        self.quality = []
        self.quality_min = 0

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
        for qval in qual:
            q = ord(qval) - 33
            self.quality.append(q)

        return self.quality

    def base_count(self):
        """-----------------------------------------------------------------------------------------
        Count the bases in the sequence that have quality > self.quality_min

        :return: dict of int, keys are bases plus 'total'
        -----------------------------------------------------------------------------------------"""
        count = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'total':0}
        for pos in range(len(self.sequence)):
            if self.quality[pos] >= self.quality_min:
                try:
                    count[self.sequence[pos]] += 1
                    count['total'] += 1
                except KeyError:
                    count[self.sequence[pos]] = 1
                    count['total'] = 1

        return count


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fastqname = '../HW1/8044.5k.fastq'
    fq = Fastq()

    fh = fq.openfile(fastqname)

    n_entry = 0
    all = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'total':0}
    while fq.next():
        n_entry += 1

        count = fq.base_count()
        for base in count:
            all[base] += count[base]

        if n_entry > 4:
            break

    # end of loop of Fastq entries
    print('{} entries read'.format(n_entry))
    print('\nAll bases')
    for base in all:
        print('{:>10s}: {}'.format(base, all[base]))

    exit(0)

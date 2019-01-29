"""=================================================================================================
Homework 2

read a fastq file and report the following

Number of sequence entries
number of bases read
number of bases with quality >= 20
ACGT count of the bases with Quality >= 20
number of bases when sequences are truncated at the first base with quality < 20
average number of bases per entry (count paired reads as separate entries).

Michael Gribskov     21 January 2019
================================================================================================="""


def quality(qstring):
    """---------------------------------------------------------------------------------------------
    Convert the quality string to a list of quality values
    quality = ord(q) - 33    where q is the character in the quality string
    the quality is -10 * log10( P(error) )

    :param qstring: character encoded quality string from fastq
    :return:
    ---------------------------------------------------------------------------------------------"""
    quality = []
    for q in qstring:
        quality.append(ord(q) - 33)

    return quality


def composition(sequence):
    """---------------------------------------------------------------------------------------------
    Return a dictionary with the count of every character in sequence

    :param sequence: DNA sequence string, may contain ambiguous bases
    :return: dictionary, indexed by bases
    ---------------------------------------------------------------------------------------------"""
    count = {}
    for base in sequence:
        if base in count:
            count[base] += 1
        else:
            count[base] = 1

    return count


def trim(sequence, qstring, qthreshold):
    """---------------------------------------------------------------------------------------------
    Truncates the sequence at the first base

    :param sequence: string, DNA sequence
    :param qstring: string, fastq quality string
    :param qthreshold: int, lowest allowed quality
    :return: string, truncated sequence string
    ---------------------------------------------------------------------------------------------"""
    q = quality(qstring)

    endpos = 0
    while q[endpos] >= qthreshold:
        endpos += 1
        if endpos == len(q):
            break

    return sequence[:endpos]


def count_hq(sequence, qstring, qthreshold):
    """---------------------------------------------------------------------------------------------

    :param sequence: string
    :param qstring: string
    :param qthreshold: int
    :return: dict
    ---------------------------------------------------------------------------------------------"""
    q = quality(qstring)

    comp = {'total': 0}
    for i in range(len(q)):
        if q[i] >= qthreshold:
            comp['total'] += 1
            if sequence[i] in comp:
                comp[sequence[i]] += 1
            else:
                comp[sequence[i]] = 1

    return comp


# --------------------------------------------------------------------------------------------------
# main program
# the if test isn't strictly necessary, but some indication that the methods are done and the
# main is beginning is very useful
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    fastqname = '../HW1/8044.5k.fastq'
    quality_threshold = 20  # pre-adaptation for future command line use

    # open the file for reading, terminate with unsuccessful status is file cannot be read
    fastq = None
    try:
        fastq = open(fastqname, 'r')
    except (IOError, OSError):
        print('Error opening file {}'.format(fastqname))
        exit(1)

    # read each line and count lines, entries, and bases
    n_line = 0
    n_entry = 0
    n_base = 0
    n_base_hq = 0
    count_base = {}
    count_base_hq = {}
    comp_base_hq = {}

    sequence = ''
    for line in fastq:

        if n_line % 4 == 1:
            # a sequence line
            n_entry += 1
            sequence = line.rstrip()
            n_base += len(sequence)

            bases = composition(sequence)
            for base in bases:
                if base in count_base:
                    count_base[base] += bases[base]
                else:
                    count_base[base] = bases[base]

        elif n_line % 4 == 3:
            # quality line
            q = quality(line.rstrip())
            trimmed = trim(sequence, line.rstrip(), quality_threshold)
            n_base_hq += len(trimmed)

            bases = composition(trimmed)
            for base in bases:
                if base in count_base_hq:
                    count_base_hq[base] += bases[base]
                else:
                    count_base_hq[base] = bases[base]

            # alternate method to count
            basecomp = count_hq(sequence, line.rstrip(), quality_threshold)
            for base in basecomp:
                if base in comp_base_hq:
                    comp_base_hq[base] += basecomp[base]
                else:
                    comp_base_hq[base] = basecomp[base]

        n_line += 1

    # end of loop over lines of file

    # report
    print('Fastq file: {}'.format(fastqname))
    print('Entries read: {}'.format(n_entry))

    print('\nTotal bases read: {}'.format(n_base))
    for base in count_base:
        print('\t{}: {}'.format(base, count_base[base]))

    print('\nHigh quality bases (Q>={}): {}'.format(quality_threshold, n_base_hq))
    for base in count_base_hq:
        print('\t{}: {}'.format(base, count_base_hq[base]))

    print('\tAverage trimmed length: {}'.format(n_base_hq / n_entry))

    print('\nHigh quality bases (Q>={}) - method 2: {}'.format(quality_threshold, n_base_hq))
    for base in comp_base_hq:
        print('\t{}: {}'.format(base, comp_base_hq[base]))

# exit with successful status
exit(0)

"""=================================================================================================
fastqCount

Counts the number of bases in a fastq formatted sequence file and reports
    total number of bases (total and by base)
    number of high quality bases with quality >= threshold (total and by base)
    number of trimmed bases (total and by base)
Trimmed bases are the bases before the first low quality base in each sequence

Michael Gribskov     30 January 2019
================================================================================================="""

def get_quality_list(qual):
    """---------------------------------------------------------------------------------------------
    Convert the quality string to a list of integer quality values.

    The ASCII characters in the quality line represent the probability that the base call at the
    corresponding position in the sequenceis incorrect. The are obtained as follows

    Q = int( -10 * log10( P(error) ) )

    Character = chr( 33 + Q )

    Q is referred to as the Phred quality, or more often as simply quality . The character '5',
    which is ASCII 53, represents quality = 20, which means an error probability of 10-2 or 0.01.
    A conventional choice for the minimum permissible quality is 20.

    Usage
        quality = get_quality_list(qual)

    :param qual: string
    :return: list of int
    ---------------------------------------------------------------------------------------------"""
    quality = []
    for qval in qual:
        q = ord(qval) - 33
        quality.append(q)

    return quality

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    fastqname = '../HW1/8044.5k.fastq'
    quality_threshold = 20  # pre-adaptation for future command line use

    # open the fastq file and read into a list of lines. terminate with unsuccessful status if
    # file cannot be read
    # Defining fq as None prevents the warning about possibly undefined variable
    fq = None
    try:
        fq = open(fastqname, 'r')
    except (IOError, OSError):
        print('Error opening file {}'.format(fastqname))
        exit(1)

    data = fq.readlines()

    # collect data for all bases, hq bases, and trimmed bases.  Use sum to store the total counts
    count_hq = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'sum': 0}
    count_all = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'sum': 0}
    count_trim = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'sum': 0}

    n_line = 0
    seq = ''
    qual = ''
    for line in data:

        if n_line % 4 == 1:
            # Save the line if it is a line a sequence line
            seq = line.rstrip()

        if n_line % 4 == 3:
            # Save the line if it is a quality line, after converting to a numeric string
            qual = line.rstrip()
            quality = get_quality_list(qual)

            # process this entry
            count_all['sum'] += len(seq)

            trimmed = False
            for pos in range(len(seq)):
                # for each sequence position, count bases
                count_all[seq[pos]] += 1

                if quality[pos] >= quality_threshold:
                    # if high quality, count HQ bases
                    count_hq[seq[pos]] += 1
                    count_hq['sum'] += 1
                    if not trimmed:
                        # do not count if a low quality base has been seen
                        count_trim[seq[pos]] += 1
                        count_trim['sum'] += 1

                else:
                    # quality is < threshold, all position after this are trimmed
                    trimmed = True

        n_line += 1

    # end of loop over data lines

    print('\nAll bases:')
    for base in count_all:
        print('\t{}: {}'.format(base, count_all[base]))

    print('\nHQ bases:')
    for base in count_hq:
        print('\t{}: {}'.format(base, count_hq[base]))

    print('\nTrimmed bases:')
    for base in count_trim:
        print('\t{}: {}'.format(base, count_trim[base]))

    exit(0)

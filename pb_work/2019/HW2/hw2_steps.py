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

    CURRENTLY  A DUMMY FUNCTION THAT RETURNS A LIST COUNTING DOWN FROM 100

    :param qual: string
    :return: list of int
    ---------------------------------------------------------------------------------------------"""
    q = 100
    quality = []
    for pos in qual:
        quality.append(q)
        q -= 1

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

            for pos in range(len(seq)):
                # for each sequence position, count bases
                count_all[seq[pos]] += 1

                if quality[pos] >= quality_threshold:
                    # if high quality, count HQ bases
                    count_hq[seq[pos]] += 1
                    count_hq['sum'] += 1

                # if not trimmed, count trimmed bases
        n_line += 1

    exit(0)

"""=================================================================================================
fastqCount

Counts the number of bases in a fastq formatted sequence file and reports
    total number of bases (total and by base)
    number of high quality bases with quality >= threshold (total and by base)
    number of trimmed bases (total and by base)
Trimmed bases are the bases before the first low quality base in each sequence

Michael Gribskov     30 January 2019
================================================================================================="""

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

    # For each line of data

        # Save the line if it is a line a sequence line
        # Save the line if it is a quality line, after converting to a numeric string

            # process this entry

                # for each sequence position
                    # count bases
                    # if high quality, count HQ bases
                    # if not trimmed, count trimmed bases

    exit(0)

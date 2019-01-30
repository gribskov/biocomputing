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

    # open the fastq file and read into a list of lines

    # For each line of data

        # Save the line if it is a line a sequence line
        # Save the line if it is a quality line, after converting to a numeric string

            # process this entry

                # for each sequence position
                    # count bases
                    # if high quality, count HQ bases
                    # if not trimmed, count trimmed bases

    exit(0)
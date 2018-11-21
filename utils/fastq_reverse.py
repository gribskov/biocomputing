"""-------------------------------------------------------------------------------------------------
Reverse complement the read sequence and quality.  this should chage "outies" to "innies"

20 November 2018     Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys


def read_fastq(fh):
    """---------------------------------------------------------------------------------------------
    Read the next 4 lines, a fastq read from the filehandle, fh.
    :param fh: open filehandle
    :return: fq, dict, keys = /title seq sep qual/
    ---------------------------------------------------------------------------------------------"""

    fq = {'title': fh.readline(), 'seq': fh.readline(), 'sep': fh.readline(), 'qual': fh.readline()}
    # fq['title'] = fh.readline()
    # fq['seq'] = fh.readline()
    # fq['sep'] = fh.readline()
    # fq['qual'] = fh.readline()

    return fq


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    f1 = sys.argv[1]
    fout = f1 + '.reverse'

    sys.stdout.write('fastq_reverse\n')
    sys.stdout.write('     fastq input: {}\n'.format(f1))
    sys.stdout.write('     fastq output: {}\n'.format(fout))

    # open input files]
    r1 = None
    try:
        r1 = open(f1, 'r')
    except IOError:
        sys.stderr.write('Unable to open fastq 1 ({})\n'.format(f1))

    # open output files
    o1 = None
    try:
        o1 = open(fout, 'w')
    except IOError:
        sys.stderr.write('Unable to open output fastq 1 ({})\n'.format(fout))

    nread = 0
    nwritten = 0
    while True:
        fastq = read_fastq(r1)
        if fastq['title']:
            break

        nread = nread + 1
        # if nread > 900:
        #     nread -= 1
        #     break
        if not nread % 1000000:
            sys.stderr.write('\n{}\t{}\n'.format(nread, nwritten))

    sys.stdout.write('\n{} reads read\n'.format(nread))
    sys.stdout.write('{} reads written to {}\n'.format(nwritten, fout))

exit(0)

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

    sys.stdout.write('fastq_trimlength\n')
    sys.stdout.write('     length: {}\n'.format(L))
    sys.stdout.write('     fastq 1: {}\n'.format(f1))
    sys.stdout.write('     fastq 2: {}\n'.format(f2))

    # open input files]
    r1 = None
    try:
        r1 = open(f1, 'r')
    except IOError:
        sys.stderr.write('Unable to open fastq 1 ({})\n'.format(f1))

    # open output files
    out1 = f1 + '.clipped'
    o1 = None
    try:
        o1 = open(out1, 'w')
    except IOError:
        sys.stderr.write('Unable to open output fastq 1 ({})\n'.format(out1))

    nread = 0
    nwritten = 0
    ndropped = 0
    limit = L + 1  # sequence strings have a newline so add one
    while True:
        fq1 = read_fastq(r1)
        if fq1['title'] == '' or fq2['title'] == '':
            break

        nread = nread + 1
        # if nread > 900:
        #     nread -= 1
        #     break
        if not nread % 1000000:
            sys.stderr.write('\n{}\t{}\t{}'.format(nread, nwritten, ndropped))

        if (len(fq1['seq']) < limit) or (len(fq2['seq']) < limit):
            ndropped += 1
            # print(len(fq1['seq']), fq1['seq'])
            # print(len(fq2['seq']), fq2['seq'])
        else:
            fq1['seq'] = fq1['seq'][:L]
            fq2['seq'] = fq2['seq'][:L]
            fq1['qual'] = fq1['qual'][:L]
            fq2['qual'] = fq2['qual'][:L]

            nwritten += 1
            o1.write('{}\n{}\n{}\n{}\n'.format(fq1['title'].rstrip(), fq1['seq'].rstrip(),
                                               fq1['sep'].rstrip(), fq1['qual'].rstrip()))
            o2.write('{}\n{}\n{}\n{}\n'.format(fq2['title'].rstrip(), fq2['seq'].rstrip(),
                                               fq2['sep'].rstrip(), fq2['qual'].rstrip()))

    sys.stdout.write('\n{} reads read\n'.format(nread))
    sys.stdout.write('{} reads written to {} and {}\n'.format(nwritten, out1, out2))
    sys.stdout.write('{} reads dropped\n'.format(ndropped))

exit(0)

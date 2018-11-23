"""-------------------------------------------------------------------------------------------------
Reverse complement the read sequence and quality.  this should chage "outies" to "innies"

20 November 2018     Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys

complement = str.maketrans('acgtunACGTUN', 'tgcaanTGCAAN')


def fastq_read(fh):
    """---------------------------------------------------------------------------------------------
    Read the next 4 lines, a fastq read from the filehandle, fh.
    :param fh: open filehandle
    :return: fq, dict, keys = /title seq sep qual/
    ---------------------------------------------------------------------------------------------"""

    fq = {'title': fh.readline().rstrip(),
          'seq': fh.readline().rstrip(),
          'sep': fh.readline().rstrip(),
          'qual': fh.readline().rstrip()}
    # fq['title'] = fh.readline()
    # fq['seq'] = fh.readline()
    # fq['sep'] = fh.readline()
    # fq['qual'] = fh.readline()

    return fq


def fastq_write(fastq, fh):
    """---------------------------------------------------------------------------------------------

    :param fastq: dictionary with fastq information: keys - title, seq, sep, qual
    :param fh: output filehandle
    :return:
    ---------------------------------------------------------------------------------------------"""
    fh.write('{}\n'.format(fastq['title']))
    fh.write('{}\n'.format(fastq['seq']))
    fh.write('{}\n'.format(fastq['sep']))
    fh.write('{}\n'.format(fastq['qual']))

    return


def sequence_revcomp(seq):
    """--------------------------------------------------------------------------------------------
    reverse complement the sequence
    :param seq:
    :return: string, reverse complement of sequence--
    ---------------------------------------------------------------------------------------------"""
    return seq.translate(complement)[::-1]


def quality_rev(qual):
    """---------------------------------------------------------------------------------------------
    Reverse the quality vector to correspond to reverse complement of sequence
    :param qual: list of ascii quality scores
    :return: list (reversed)
    ---------------------------------------------------------------------------------------------"""
    return qual[::-1]


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    f1 = sys.argv[1]
    fout = f1 + '.reverse'

    sys.stdout.write('fastq_reverse\n')
    sys.stdout.write('     fastq input: {}\n'.format(f1))
    sys.stdout.write('     fastq output: {}\n'.format(fout))

    # open input file
    r1 = None
    try:
        r1 = open(f1, 'r')
    except IOError:
        sys.stderr.write('Unable to open fastq 1 ({})\n'.format(f1))

    # open output file
    o1 = None
    try:
        o1 = open(fout, 'w')
    except IOError:
        sys.stderr.write('Unable to open output fastq 1 ({})\n'.format(fout))

    nread = 0
    nwritten = 0
    while True:
        fastq = fastq_read(r1)
        if not fastq['title']:
            break

        nwritten += 1
        fastq['seq'] = sequence_revcomp(fastq['seq'])
        fastq['qual'] = quality_rev(fastq['qual'])
        fastq_write(fastq, o1)

        nread = nread + 1
        if nread > 10:
            break

        if not nread % 1000000:
            sys.stderr.write('\n{}\t{}\n'.format(nread, nwritten))

    sys.stdout.write('\n{} reads read\n'.format(nread))
    sys.stdout.write('{} reads written to {}\n'.format(nwritten, fout))

exit(0)

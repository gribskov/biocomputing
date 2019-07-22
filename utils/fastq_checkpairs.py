"""=================================================================================================
fastq_checkpairs

check that the two fastq files are properly paired

13 march 2019   Michael Gribskov
================================================================================================="""
import sys
import re
from sequence.fastq import Fastq


def sra_id(fastq):
    """---------------------------------------------------------------------------------------------

    :param fastq:
    :return:
    ---------------------------------------------------------------------------------------------"""
    sranum = re.compile(r'[^\.]+\.(\d+)')

    seqnum = sranum.search(fastq.id)
    # print('Match found: ', seqnum.group(1))
    num = 0
    try:
        num = seqnum.group(1)
    except:
        sys.stdout.write('error\n')
        fastq.write(sys.stdout)
        exit(2)

    return num


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    SRA = True  # for sequences from NCBI short read archive

    f1 = sys.argv[1]
    f2 = sys.argv[2]
    fq1 = Fastq(filename=f1)
    fq2 = Fastq(filename=f2)

    f1out = f1 + '.fixed'
    f2out = f2 + '.fixed'
    try:
        out1 = open(f1out, 'w')
        out2 = open(f2out, 'w')
    except (IOError,OSError):
        sys.stderr.write( 'unable to open output files({}, {})'.format(f1out, f2out))
        exit(1)

    sys.stdout.write('fastq_checkpairs\n')
    sys.stdout.write('     fastq input 1: {}\n'.format(f1))
    sys.stdout.write('     fastq input 2: {}\n'.format(f2))

    sys.stdout.write('     fixed fastq output: {}\n'.format(f1out))
    sys.stdout.write('     fixed fastq output: {}\n'.format(f2out))

    OK = True
    n_read = 0
    while OK:
        fq1.next()
        fq2.next()

        if not fq1.quality:
            sys.stderr.write('No entry read from {}'.format(f1))
            exit(2)
        if not fq2.quality:
            sys.stderr.write('No entry read from {}'.format(f2))
            exit(2)

        if not fq1.check:
            sys.stdout.write('{} error: {}'.format(fq1, fq1.error))
            OK = False
        if not fq1.check:
            sys.stdout.write('{} error: {}'.format(fq2, fq2.error))
            OK = False

        if not OK:
            break

        n_read += 1
        print("{}".format(n_read))

        # for SRA files, check to see the read number matches
        match = False
        if SRA:
            if sra_id(fq1) == sra_id(fq2):
                match = True
        if match:
            fq1.write(out1)
            fq2.write(out2)

exit(0)

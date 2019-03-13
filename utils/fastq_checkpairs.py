"""=================================================================================================
fastq_checkpairs

check that the two fastq files are properly paired

13 march 2019   Michael Gribskov
================================================================================================="""
import sys
from sequence.fastq import Fastq

# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    f1 = sys.argv[1]
    f2 = sys.argv[2]
    fq1 = Fastq(filename=f1)
    fq2 = Fastq(filename=f2)

    f1out = f1 + '.fixed'
    f2out = f2 + '.fixed'

    sys.stdout.write('fastq_checkpairs\n')
    sys.stdout.write('     fastq input 1: {}\n'.format(f1))
    sys.stdout.write('     fastq input 2: {}\n'.format(f2))

    sys.stdout.write('     fixed fastq output: {}\n'.format(f1out))
    sys.stdout.write('     fixed fastq output: {}\n'.format(f2out))

    OK = True
    n_read = 0
    while OK:
        fq1.read()
        fq2.read()

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

        sys.stdout.write(fq1.write())
        sys.stdout.write(f12.write())

exit(0)
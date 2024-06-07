"""=================================================================================================
Convert a fasta file to a fastq file with dummy quality
assumes the sequence is on a single line

Michael Gribskov     07 June 2024
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    qual = sys.argv[2]
    fa_name = sys.argv[1]
    fa = open(fa_name, 'r')

    fq_name = fa_name.replace('.fasta', '.fq')
    fq_name = fa_name.replace('.fa', '.fq')
    fq = open( fq_name, 'w')

    fa_n = 0
    fq_n = 0
    for line in fa:

        if line.startswith('>'):
            fq.write(line.rstrip().replace('>','@', 1))
            fa_n += 1
        else:
            seq = line.rstrip()
            fq.write(f' length={len(seq)}\n')
            fq.write(seq)
            fq.write('\n+\n')
            fq.write(f'{qual*len(seq)}\n')
            fq_n += 1

    print(f'{fa_n} fasta entries read')
    print(f'{fq_n} fastq entries written')

    exit(0)
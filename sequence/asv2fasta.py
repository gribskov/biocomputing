"""=================================================================================================
convert a list of ASVs and sequences from Qiime to a fasta file for blast searching

Michael Gribskov     15 April 2024
================================================================================================="""
import sys
# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    asvfile = sys.argv[1]
    asv = open(asvfile, 'r')

    fastafile = sys.argv[2]
    fasta = open(fastafile, 'w')

    nseq = 0
    for sequence in asv:
        if not sequence:
            continue
        nseq += 1
        asv, seq = sequence.rstrip().split()
        fasta.write(f'{asv}\n{seq}\n')

    asv.close()
    fasta.close()

    print(f'{nseq} reformatted')

    exit(0)

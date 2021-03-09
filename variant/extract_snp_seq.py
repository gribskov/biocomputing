"""=================================================================================================
Extract SNP sequences from the genome given the chromosome and base position
+/- pad bases are extracted around the position


Michael Gribskov     08 March 2021
================================================================================================="""
from sequence.fasta import Fasta


def read_snps_tabular(filename):
    """---------------------------------------------------------------------------------------------
    format is
    ID<tab>Chromosome number<tab>base

    :param filename: string, file with SNPs
    :return: list of dict with keys id, chr, pos, in_25k_subset
    ---------------------------------------------------------------------------------------------"""
    try:
        snp = open(filename, 'r')
    except (IOError, OSError):
        print('unable to open SNP file ({})'.format(filename))
        exit(1)

    snplist = {}
    maxpos = {}
    for line in snp:
        if line.startswith('SNP '):
            # skip header
            print('skipping line {}'.format(line))
            continue

        field = line.rstrip().split('\t')

        # check if in 25K list
        shortlist = False
        try:
            if field[3] == 'yes':
                shortlist = True
        except IndexError:
            pass

        if field[1] == '18':
            # in the genome unplaced scaffolds are assigned chromosome 00
            field[1] = '00'
        chr = 'Chr{:02d}'.format(int(field[1]))

        if chr in snplist:
            snplist[chr].append({'id': field[0], 'pos': int(field[2]), 'in_25k': shortlist})
            maxpos[chr] = max( maxpos[chr], int(field[2]))
        else:
            snplist[chr] = [{'id': field[0], 'pos': int(field[2]), 'in_25k': shortlist}]
            maxpos[chr] = int(field[2])

    snp.close()

    return snplist, maxpos


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read SNPS
    robust_snp = 'C:/Users/michael/Desktop/apple/235000_robust_snps.txt'
    snplist, maxpos = read_snps_tabular(robust_snp)
    for chr in snplist:
        print('{} {} snps max: {}'.format(chr, len(snplist[chr]), maxpos[chr]))

    # read genome and match, one sequence at a time
    fastafile = 'C:/Users/michael/Desktop/apple/GDDH13_1-1_formatted.fasta'
    fasta = Fasta()
    fasta.open(fastafile)

    bases = 0
    sequence = {}
    for line in fasta.fh:
        if line.startswith('>'):
            try:
                id, doc = line.rstrip().split(' ')
            except ValueError:
                id = line.rstrip()

            id = id.replace('>', '')
            sequence[id] = 0

        else:
            sequence[id] += len(line.rstrip())

    for chr in sorted(sequence):
        print('{}\t{} bases\t{} snps\t{} max'.format(chr,sequence[chr], len(snplist[chr]), maxpos[chr]))

    exit(0)

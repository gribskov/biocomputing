"""=================================================================================================
Extract SNP sequences from the genome given the chromosome and base position
+/- pad bases are extracted around the position


Michael Gribskov     08 March 2021
================================================================================================="""
import sys
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
            maxpos[chr] = max(maxpos[chr], int(field[2]))
        else:
            snplist[chr] = [{'id': field[0], 'pos': int(field[2]), 'in_25k': shortlist}]
            maxpos[chr] = int(field[2])

    snp.close()

    return snplist, maxpos


def read_snps_gff(filename):
    """---------------------------------------------------------------------------------------------
    format is
    Chr11   IRHS2017        feature 34037208        34037208        1.0     +       .       ID=FB_AFFY_7878392;ax=AX-115460897;affx=Affx-114118547
    Chr15   IRHS2017        feature 50582376        50582376        1.0     -       .       ID=FB_AFFY_7878301;ax=AX-115460894;affx=Affx-113780418

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
        field = line.rstrip().split('\t')

        # get SNP ID from attributes
        token = field[8].split(';')
        for pairs in token:
            if pairs.startswith('ID='):
                tag, affyid = pairs.split('=')

        chr = field[0]
        if chr in snplist:
            snplist[chr].append({'id': affyid, 'pos': int(field[3]), 'strand': field[6]})
            maxpos[chr] = max(maxpos[chr], int(field[3]))
        else:
            snplist[chr] = [{'id': affyid, 'pos': int(field[3]), 'strand': field[6]}]
            maxpos[chr] = int(field[3])

    for chr in snplist:
        #sort by position in each chromosome
        snplist[chr].sort(key=lambda s:s['pos'])

    snp.close()

    return snplist, maxpos


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read SNPS
    # this version readstabular format
    snp_file = 'C:/Users/michael/Desktop/apple/235000_robust_snps.txt'
    snplist, maxpos = read_snps_tabular(snp_file)

    # this version reads gff format
    # snp_file = 'C:/Users/michael/Desktop/apple/GDDH13_1-1_SNPs.gff3'
    # snplist, maxpos = read_snps_gff(snp_file)

    for chr in snplist:
        print('{} {} snps max: {}'.format(chr, len(snplist[chr]), maxpos[chr]))
    print()

    # output file from command line
    out = open(sys.argv[1], 'w')

    # read genome and match, one sequence at a time
    fastafile = 'C:/Users/michael/Desktop/apple/GDDH13_1-1_formatted.fasta'
    fasta = Fasta()
    fasta.open(fastafile)

    pad = 25
    window = 2 * pad + 1

    bases = 0
    seqlen = {}
    sequence = ''
    seqbegin = 0
    seqend = 0
    snpcount = 0
    for line in fasta.fh:
        line = line.rstrip()
        if line.startswith('>'):
            try:
                id, doc = line.split(' ')
            except ValueError:
                id = line

            id = id.replace('>', '')
            seqlen[id] = 0
            seqbegin = 0
            seqend = 0
            snpset = (s for s in snplist[id])
            try:
                snp = next(snpset)
                fastforward = False
            except StopIteration:
                break
            snp_pos = snp['pos']
            snp_id = snp['id']

        else:
            if fastforward:
                continue

            l = len(line)
            seqlen[id] += l
            sequence += line
            seqend = seqbegin + len(sequence)

            if snp_pos > seqend + pad:
                # snp pos is not in this sequence
                seqbegin += len(sequence)
                sequence = ''
                continue

            while snp_pos + pad < seqend:
                start = snp_pos - seqbegin - pad - 1
                stop = start + window
                snpseq = sequence[start:stop]
                # out.write('{}\t{}\t{}\t{}\n'.format(id, snp_pos, snpseq, snp_id))
                out.write('>{} {}:{}\n{}\n'.format(snp_id, id, snp_pos, snpseq))
                snpcount += 1
                try:
                    snp = next(snpset)
                    snp_pos = snp['pos']
                    snp_id = snp['id']
                except StopIteration:
                    next_chr = True
                    fastforward = True
                    break

    print('{} SNPs written'.format(snpcount))

    exit(0)

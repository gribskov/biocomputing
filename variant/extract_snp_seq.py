"""=================================================================================================
Extract SNP sequences from the genome given the chromosome and base position
+/- pad bases are extracted around the position


Michael Gribskov     08 March 2021
================================================================================================="""


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

    snplist = []
    for line in snp:
        if line.startswith('SNP'):
            # skip header
            continue

        field = line.rstrip().split('\t')

        # check if in 25K list
        shortlist = False
        try:
            if field[3] == 'yes':
                shortlist = True
        except IndexError:
            pass

        snplist.append({'id': field[0], 'chr': field[1], 'pos': int(field[2]), 'in_25k':shortlist})

    snp.close()
    print('{}SNPs read from {}'.format(len(snplist), filename))

    return snplist


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read SNPS
    robust_snp = 'C:/Users/michael/Desktop/apple/235000_robust_snps.txt'
    read_snps_tabular(robust_snp)


    # read genome and match, one sequence at a time
    exit(0)

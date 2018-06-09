"""=================================================================================================
Merge multiple count files from RSEM into a single file.
input (tab delimited):
gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
TRINITY_DN10023 TRINITY_DN10023_c0_g1_i1        522.00  422.34  0.00    0.00    0.00
TRINITY_DN10027 TRINITY_DN10027_c0_g1_i1        552.00  452.34  0.00    0.00    0.00
TRINITY_DN1003  TRINITY_DN1003_c0_g1_i1,TRINITY_DN1003_c0_g2_i1 335.00  235.34  0.00    0.00    0.00
TRINITY_DN10038 TRINITY_DN10038_c0_g1_i1,TRINITY_DN10038_c0_g1_i2       237.50  137.84  0.00    0.00    0.00

This versiuon is tailored for names that look like 024236_TE1.rsem.genes.results where the second
token, TE1, indicates replicate 1 of sample TE

usage:
    rsem_merge_counts.py <rsem_result>

Michael Gribskov     09 June 2018
================================================================================================="""
import sys
import os
import glob


def select_by_group(value, group, minval):
    """---------------------------------------------------------------------------------------------
    returns true if any group has equal or more than minval counts
    :param value: dict of count values
    :param group: list of data groups
    :param minval: integer, minimum value must be greater for any group
    :return:
    ---------------------------------------------------------------------------------------------"""
    if value['_total'] < minval:
        return False

    for g in group:
        if value[g] >= minval:
            return True

    return False


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    junk = []
    data = {}
    group = []  # names of sample groups
    lastgene = ''
    minval = int(sys.argv[2])

    for wild in glob.iglob(sys.argv[1]):
        path, file = os.path.split(wild)
        sys.stderr.write('processing file: {}\n'.format(file))

        # get sample
        sample, junk = file.split('.', maxsplit=1)
        samplenum, samplename = sample.split('_')
        samplegroup = samplename.rstrip('0123456789')
        sys.stderr.write('    sample number: {}    name: {}    group: {}\n'.format(
            samplenum, samplename, samplegroup))

        if samplename in data:
            sys.stderr.write('Error: duplicate sample name({})\n'.format(samplename))
            sys.stderr.write('    Skipping sample\n')
            continue

        try:
            rsem = open(wild, 'r')
        except:
            sys.stderr.write('Error opening rsem results file ({})'.format(wild))
            continue

        # save in list of sample groups if new
        if samplegroup not in group:
            group.append(samplegroup)

        nline = 0
        nvalue = 0
        nzero = 0
        count_sum = 0.0
        icount_sum = 0

        # first line is header, skip
        line = rsem.readline()

        # column for each sample and for sum over sample group
        for line in rsem:
            nline += 1
            # fields are gene_id, transcript_id(s), length, effective_length, expected_count, TPM, FPKM
            field = line.split()
            gene = field[0]
            count = float(field[4])
            icount = round(count)

            count_sum += count
            icount_sum += icount

            if count == 0:
                nzero += 1
            else:
                nvalue += 1

            if gene not in data:
                data[gene] = {}
                lastgene = gene

            data[gene][samplename] = icount
            for set in [samplegroup, '_total']:
                if set in data[gene]:
                    data[gene][set] += icount
                else:
                    data[gene][set] = icount

        sys.stderr.write('    lines: {}\n'.format(nline))
        sys.stderr.write('    values: {}\n'.format(nvalue))
        sys.stderr.write('    zeroes: {}\n'.format(nzero))
        sys.stderr.write('    sum(float): {:.2f}\n'.format(count_sum))
        sys.stderr.write('    sum(int): {}\n'.format(icount_sum))

    # end of loop over files

    # tabular result
    sys.stdout.write('#gene')
    for gene in sorted(data[lastgene].keys()):
        sys.stdout.write('\t{}'.format(gene))
    sys.stdout.write('\n')

    nselect = 0
    for gene in sorted(data.keys(), key=lambda x: data[x]['_total']):
        if select_by_group(data[gene], group, minval):
            sys.stdout.write('{}'.format(gene))
            for set in sorted(data[gene].keys()):
                sys.stdout.write('\t{}'.format(data[gene][set]))
            sys.stdout.write('\n')
            nselect += 1

    sys.stderr.write('\n{} genes selected\n'.format(nselect))

    exit(0)

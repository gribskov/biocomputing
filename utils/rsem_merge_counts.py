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
import glob

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    junk = []
    data = {} = {}
    for file in glob.iglob(sys.argv[1])
        sys.stderr.write('processing file: {}\n'.format(file))

        # get sample
        sample, junk = file.split('.')
        samplenum, samplename = sample.split('_')
        samplegroup = samplename.rstrip('0123456789')
        sys.stderr.write('    sample number: {}    name: {}    group: {}\n'.format(
            samplenum, samplename, samplegroup))

        if samplename in data:
            sys.stderr.write('Error: duplicate sample name({})\n'.format(samplename))
            sys.stderr.write('    Skipping sample\n')
            continue

        try:
            rsem = open(file, 'r')
        except:
            sys.stderr.write('Error opening rsem results file ({})'.format(trinity_file))

        nline = 0
        nvalue = 0
        nzero = 0
        count_sum = 0.0
        icount_sum = 0

        # first line is header, skip
        line = rsem.readline()

        # column for each sample and for sum over sample group
        data[samplename] = {}
        if samplegroup not in data:
            data[group] = {}
            
        for line in rsem:
            nline += 1
            # fields are gene_id, transcript_id(s), length, effective_length, expected_count, TPM, FPKM
            field = line.split()
            gene = field[0]
            count = field[4]
            icount = round(count)
            
            count_sum += count
            icount_sum += icount
            
            if count == 0:
                nzero += 1
            else:
                nvalue += 1
                
            data[gene][samplename] = count
            for set in [samplegroup, 'total']:
                if set in data[gene]:
                    data[gene][set] += count
                else:
                    data[gene][set] = count

        sys.stderr.write('    lines: {}\n'.format(nline))
        sys.stderr.write('    values: {}\n'.format(nvalue))
        sys.stderr.write('    zeroes: {}\n'.format(nzero))
        sys.stderr.write('    sum(float): {}\n'.format(count_sum))
        sys.stderr.write('    sum(int): {}\n'.format(icount_sum))
            
    # end of loop over files

    # tabular result
    for gene in sorted(data.keys()):
        sys.stdout.write(gene)
        for set in sorted[data[gene]]:
            sys.stdout.write('\t{}'.format(data[gene][set])
        sys.stdout.write('\n')

    exit(0)
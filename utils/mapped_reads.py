"""=================================================================================================
Count the number of mapped reads

Michael Gribskov     08 July 2018
================================================================================================="""
import sys


def trinity_shortname(name, level='c'):
    nfield = {
        'cluster': 2,
        'l': 2,
        'component': 3,
        'c': 3,
        'group': 4,
        'gene': 4,
        'g': 4,
        'isoform': 5,
        'i': 5
    }
    field = name.split('_')

    return '_'.join(field[:nfield[level]])


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    count_file = sys.argv[1]
    level = sys.argv[2]

    sys.stderr.write('mapped_reads\n\tinput: {}\n\tlevel: {}\n\n'.format(count_file, level))

    try:
        count = open(count_file, 'r')
    except:
        sys.stderr.write('Error opening count file ({})'.format(count_file))
        exit(1)

    header = count.readline()
    sys.stdout.write(header)

    nline = 0
    nmapped = 0
    ngene = 0
    gene = {}
    for line in count:
        nline += 1

        name, l, eff_l, tpm, c = line.rstrip().split()
        l = int(l)
        eff_l = float(eff_l)
        tpm = float(tpm)
        c = float(c)
        nmapped += c
        short = trinity_shortname(name, level)

        print('{}\t{}'.format(short, name))
        if short in gene:
            gene[short]['tpm'] += tpm
            gene[short]['count'] += c
            gene[short]['length'] += l
            gene[short]['effective_length'] += eff_l
            gene[short]['n'] += 1

        else:
            ngene += 1
            sys.stderr.write('{}\t{}\n'.format(ngene, short))
            gene[short] = {}
            gene[short]['tpm'] = tpm
            gene[short]['count'] = c
            gene[short]['length'] = l
            gene[short]['effective_length'] = eff_l
            gene[short]['n'] = 1

    sys.stderr.write('\n{:.2f} read counts read\n{} isoforms\n{} genes\n'.format(nmapped, nline, ngene))

    #  write out results
    for g in sorted(gene):
        l_ave = gene[g]['length'] / gene[g]['n']
        effl_ave = gene[g]['effective_length'] / gene[g]['n']
        sys.stdout.write(
            '{}\t{:.2f}\t{:.2f}\t{:.3f}\t{:.3f}\n'.format(g, l_ave, effl_ave, gene[g]['tpm'], gene[g]['count']))

exit(0)

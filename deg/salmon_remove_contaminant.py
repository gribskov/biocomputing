####################################################################################################
# Use for trinity gene counts. Remove a list of contaminants identified by blast
#
# Michael Gribskov  1 December 2024
####################################################################################################
import sys


def contamread(infilename):
    """---------------------------------------------------------------------------------------------
    The input file is from blast_filter.py. There are two sections, the first with a list of
    trinity bundles to remove, and the second with a list of taxa. This function reads the first
    section and returns a list

    the first section looks like: (numbers are number of bundle components with blast hits
    Viral assemblies - total found = 1934
    69	DN2
    48	DN66978
    ...
     1	DN21006
     1	DN21066

    Viral taxa

    :param infilename: string    path to contaminant file
    :return: list of string      trinity bundles considered to be contaminants
    ---------------------------------------------------------------------------------------------"""
    try:
        infile = open(infilename, 'r')
    except OSError:
        print(f'Error opening contaminant file ({infilename}) on ARGV[2]')
        exit[2]

    title = infile.readline()
    blist = []
    n = 0
    for line in infile:
        if line.isspace():
            break

        count, bundle = line.rstrip().split()
        n += 1
        blist.append(bundle)

    return blist


def print_counts(outfname, argnum, colnames, data, ctmin):
    """---------------------------------------------------------------------------------------------
    
    :param outfname: string         output file name
    :param colnames:list of string  names of data columns
    :param argnum: int              field on command line and error number
    :param data: dict of list       data columns for each trinity group
    :param ctmin: int               minimum number of counts to retain
    :return: True
    ---------------------------------------------------------------------------------------------"""
    try:
        outfile = open(outfname, 'w')
    except OSError:
        print(f'Error opening salmon count file ({outfname}) on ARGV[{argnum}]')
        exit(argnum)

    outfile.write(f'gene\t{'\t'.join(f'{n}' for n in colnames)}\n')
    dropn = 0
    dropct = 0
    writtenn = 0
    writtenct = 0
    for gene in data:
        if gene == 'total':
            continue

        sum = 0
        for value in data[gene]:
            sum += value

        if sum < ctmin:
            dropn += 1
            dropct += sum
            continue

        outfile.write(f'{gene}\t{'\t'.join(f'{n:10.0f}' for n in data[gene])}\n')
        writtenct += sum
        writtenn += 1

    outfile.close()
    print(f'{writtenn + dropn} merged genes examined')
    print(f'{writtenn} merged genes with  {writtenct:10.0f} counts written to {outfname}')
    print(f'{dropn} merged genes with {dropct:10.0f} counts dropped at threshold={ctmin}')

    return True


####################################################################################################
# main
####################################################################################################
if __name__ == '__main__':
    level = int(sys.argv[1])
    remove = sys.argv[2]
    dirty = sys.argv[3]
    clean = sys.argv[4]
    contam = sys.argv[5]

    # read in contaminants
    contamlist = contamread(remove)
    print(f'{len(contamlist)} bundles IDs read from {remove}')

    # open salmon file and read column header, remove the string .salmon from column titles
    # format:
    # Name	C16T00R1.salmon	C16T00R2.salmon	C16T00R3.salmon	C16T06R1.salmon	C16T06R2.salmon	C16T06R3.salmon	C16T48R1.salmon	...
    # TRINITY_DN29267_c0_g1_i1	0	2	0	0	0	0	7	4	3	0	0	9	2	3.057	5	0	2	1	2	3	2	1	2	0
    # TRINITY_DN29246_c0_g2_i4	0	0	0	3.968	1	1	2.23	0	0	0	0	5	0	6	0	0	2	0	0	0	1	1.001	0	0
    # TRINITY_DN29246_c0_g2_i3	0	0	1	0	1	0	14.321	2	7	6	2	0	0	0	0	5	0	3	1	0	0	0	0	0
    try:
        ctfile = open(dirty, 'r')
    except OSError:
        print(f'Error opening salmon count file ({dirty}) on ARGV[3]')
        exit[3]

    colname = []
    line = ctfile.readline()
    field = line.rstrip().split('\t')
    colname = [field[x].replace('.salmon', '') for x in range(1, len(field))]
    print(f'{len(colname)} columns read from {dirty}')

    # for each row in salmon
    # sum all counts at desired level and write to either clean or contaminant output
    good = {'total': [0 for _ in range(len(colname))]}
    bad = {'total': [0 for _ in range(len(colname))]}
    genen = 0
    for line in ctfile:
        if line.isspace():
            continue

        genen += 1
        field = line.rstrip().split('\t')
        trinity = field[0].split('_')
        tlevel = '_'.join(trinity[1:level + 1])  # this is the level for gathering counts

        if trinity[1] in contamlist:
            if tlevel not in bad:
                bad[tlevel] = [0 for _ in range(len(colname))]
            for i in range(1, len(colname) + 1):
                bad[tlevel][i - 1] += float(field[i])
                bad['total'][i - 1] += float(field[i])
        else:
            if tlevel not in good:
                good[tlevel] = [0 for _ in range(len(colname))]
            for i in range(1, len(colname) + 1):
                good[tlevel][i - 1] += float(field[i])
                good['total'][i - 1] += float(field[i])

    print(f'{genen} trinity predicted isoforms processed\n')
    print(f'Total Counts by sample')
    print(f'\t\t{'\t'.join(f'{n}' for n in colname)}')
    print(f'good\t{'\t'.join(f'{n:10.0f}' for n in good['total'])}')
    print(f' bad\t{'\t'.join(f'{n:10.0f}' for n in bad['total'])}\n')

    # print out the clean and contaminant files with merged counts
    print(f'Cleaned and merged counts written to {clean}')
    print_counts(clean, 3, colname, good, 100)
    print(f'\nContaminant merged counts written to {contam}')
    print_counts(contam, 4, colname, bad, 100)

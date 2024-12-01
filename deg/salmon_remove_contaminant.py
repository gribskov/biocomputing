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


####################################################################################################
# main
####################################################################################################
if __name__ == '__main__':
    level = sys.argv[1]
    dirty = sys.argv[2]
    clean = sys.argv[3]
    contam = sys.argv[4]

    # read in contaminants
    contam = contamread(dirty)
    print(f'{len(contam)} bundles IDs read from {dirty}')

    # for each row in salmon
    # sum all counts at desired level and write to either clean or contaminant output

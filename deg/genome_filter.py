"""=================================================================================================
For a set of genes compared to a reference genome using blastn, separate into matchng and not
matching sets.
this can be used to removed contaminants IF YOU ARE SURE YOUR REFERENCE GENOME IS COMPLETE

Michael Gribskov     18 April 2025
================================================================================================="""
import sys
from collections import defaultdict


def fasta_ids(fastafname):
    """---------------------------------------------------------------------------------------------
    Read just the sequence IDs from a fasta file
    :param fastafname: string       path to fasta nucleotide sequence file
    :return: list                   IDs of all sequences
    ---------------------------------------------------------------------------------------------"""
    fasta = open(fastafname, 'r')
    ids = []
    for line in fasta:
        if line.startswith('>'):
            # title line
            field = line.split(' ')
            ids.append(field[0].lstrip('>'))

    fasta.close()
    return ids


def blastfmt7_ids(blastfname):
    """---------------------------------------------------------------------------------------------
    Return list of all IDs, and IDs that have no hits in the search. this relies on th
    :param blastfname:
    :return:
    ---------------------------------------------------------------------------------------------"""
    blast = open(blastfname, 'r')

    ids = []
    missing = []
    n = 0
    for line in blast:
        if line.startswith('# Query:'):
            # query sequence ID
            field = line.split(' ')
            thisid = field[2]
            ids.append(thisid)
        elif line.find('hits found') > -1:
            field = line.split(' ')
            if field[1] == '0':
                missing.append(thisid)

            # debugging
            n += 1
            if n > 10000:
                break

    return ids, missing


def trinity_filter(ids, level=3, remove_trinity=False):
    """---------------------------------------------------------------------------------------------
    return a new list with the trinity IDs trimmed to the indicated level:
    level 4: all parts, TRINITY_DN221212_c0_g1_i1
    level 3: gene, TRINITY_DN221212_c0_g1
    level 2: component, TRINITY_DN221212_c0
    level 1: bundle, TRINITY_DN221212
    The names in the new list are unique

    :param ids: list                    strings with trinity sequence IDs
    :param level: int                   level to trim trinity IDs, default is gene level
    :param remove_trinity: boolean      remove the use TRINITY_ prefix if true
    :return: list                       strings with trinity IDs trimmed at the desired level
    ---------------------------------------------------------------------------------------------"""
    idhash = defaultdict(int)
    first = 0
    last = level + 1
    if remove_trinity:
        first = 1

    for name in ids:
        field = name.split('_')
        thisid = '_'.join(field[first:last])
        idhash[thisid] += 1

    return idhash.keys()


# --------------------------------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # id_fasta = fasta_ids(sys.argv[1])
    # print(f'fasta sequences read: {len(id_fasta)}')

    ids, id_missing = blastfmt7_ids(sys.argv[2])
    print(f'ids read: {len(ids)}')
    print(f'missing ids: {len(id_missing)}')
    ids = trinity_filter(ids, level=2, remove_trinity=True)

    exit(0)

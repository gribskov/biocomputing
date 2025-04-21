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
    Return list of query IDs that have hits, and IDs that have no hits in the search. this relies on
    the Blast format=7 which lists queries with no hits

    :param blastfname:
    :return:
    ---------------------------------------------------------------------------------------------"""
    blast = open(blastfname, 'r')

    found = []
    missing = []
    thisid = None
    # n = 0
    for line in blast:
        if line.startswith('# Query:'):
            # query sequence ID
            field = line.split(' ')
            thisid = field[2]
            # found.append(thisid)
        elif line.find('hits found') > -1:
            field = line.split(' ')
            if field[1] == '0':
                missing.append(thisid)
            else:
                found.append(thisid)

        # # debugging
        # n += 1
        # if n > 10000:
        #     break

    return found, missing


def trinity_filter_id(tid, level=3, remove_trinity=False):
    """---------------------------------------------------------------------------------------------
    filter a trinity id to truncate it at the desired level, and optionally drop the redundant
    TRINITY_ prefix

    level 4: all parts, TRINITY_DN221212_c0_g1_i1
    level 3: gene,      TRINITY_DN221212_c0_g1
    level 2: component, TRINITY_DN221212_c0
    level 1: bundle,    TRINITY_DN221212

    :param tid: string                   trinity sequence ID
    :param level: int                   level to trim trinity IDs, default is gene level
    :param remove_trinity: boolean      remove the use TRINITY_ prefix if true
    :return: list                       strings with trinity IDs trimmed at the desired level
    ---------------------------------------------------------------------------------------------"""
    first = 0
    last = level + 1
    if remove_trinity:
        first = 1

    field = tid.split('_')
    return '_'.join(field[first:last])

def trinity_filter_list(found, missing, level=3, remove_trinity=False):
    """---------------------------------------------------------------------------------------------
    return new found and missing lists with the trinity IDs trimmed to the indicated level:
    level 4: all parts, TRINITY_DN221212_c0_g1_i1
    level 3: gene,      TRINITY_DN221212_c0_g1
    level 2: component, TRINITY_DN221212_c0
    level 1: bundle,    TRINITY_DN221212
    The names in the new list are unique
    some sequences may have both hits and no-hits results for different isoforms, in this case, all
    of the sequences are included in the found list

    :param found: list                  strings with trinity sequence IDs found
    :param missing: list                strings with trinity sequences IDs not found
    :param level: int                   level to trim trinity IDs, default is gene level
    :param remove_trinity: boolean      remove the use TRINITY_ prefix if true
    :return: list                       strings with trinity IDs trimmed at the desired level
    ---------------------------------------------------------------------------------------------"""
    idhash = defaultdict(bool)

    for name in missing:
        idhash[trinity_filter_id(name, level, remove_trinity)] = False

    for name in found:
        idhash[trinity_filter_id(name, level, remove_trinity)] = True

    found = [k for k in idhash if idhash[k]]
    missing = [k for k in idhash if not idhash[k]]

    return found, missing


def getfasta(fastafname, level=3, remove_trinity=True, idfilter=lambda x:x):
    """---------------------------------------------------------------------------------------------
    generator that yields the next fasta sequence in the file

    :param fastafname: string       path to fasta file
    :param level: int               level to truncate trinity IDs
    :param remove_trinity: bool     if true remove 'TRINITY_' from ID
    :param idfilter: func           function to filter trinity IDs (e.g., to remove TRINITY_)
    :yield: dict                    id, doc, sequence (sequence is not stripped)
    ---------------------------------------------------------------------------------------------"""
    fasta = open(fastafname, 'r')

    entry = {'id': '', 'doc': '', 'sequence': ''}
    for line in fasta:
        if line.startswith('>'):
            if entry['sequence']:
                yield entry
                entry = {'id': '', 'doc': '', 'sequence': ''}
            # title line
            field = line.rstrip().split(' ')
            entry['id'] = idfilter(field[0].lstrip('>'), level, remove_trinity)
            entry['doc'] = ''
            if len(field) > 1:
                entry['doc'] = field[1]
            entry['sequence'] = ''
        else:
            entry['sequence'] += line

    if entry['sequence']:
        yield entry

    fasta.close()
    return


# --------------------------------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    trinity_level = 3
    print(f'genomefilter.py')
    print(f'trinity level: {trinity_level}')
    print(f'blast search (trinity): {sys.argv[2]}')
    print(f'fasta sequence file (trinity): {sys.argv[3]}')
    print(f'blast search(genome) {sys.argv[2]}')

    found, nohits = blastfmt7_ids(sys.argv[2])
    print(f'trinity ids: {len(found) + len(nohits)}')
    print(f'trinity ids with hits: {len(found)}')
    print(f'trinity ids with no hits {len(nohits)}')

    found, missing = trinity_filter_list(found, nohits, level=trinity_level, remove_trinity=True)
    print(f'\nAfter filtering at trinity level={trinity_level}:')
    print(f'ids matched: {len(found)}')
    print(f'ids not matched: {len(missing)}')

    # write sequences to separate outputs
    basename = sys.argv[1].replace('.fasta', '').replace('.fa', '')
    foundout = f'{basename}.found.fa'
    foundfh = open(foundout, 'w')
    missingout = f'{basename}.missing.fa'
    missingfh = open(missingout, 'w')

    unknown = []
    found_n = 0
    missing_n = 0
    for fasta in getfasta(sys.argv[3], level=trinity_level, idfilter=trinity_filter_id):
        if fasta['id'] in missing:
            missingfh.write(f">{fasta['id']} {fasta['doc']}\n{fasta['sequence']}")
            missing_n += 1
        elif fasta['id'] in found:
            foundfh.write(f">{fasta['id']} {fasta['doc']}\n{fasta['sequence']}")
            found_n += 1
        else:
            unknown.append(fasta['id'])

    print(f'\n{found_n} sequences written to {foundout}')
    print(f'{missing_n} sequences written to {missingout}')
    if unknown:
        print(f'{len(unknown)} unclassified sequences')

    foundfh.close()
    missingfh.close()

    exit(0)

"""=================================================================================================
filter a trinity assembly using a diamond/blastx search

Michael Gribskov     01 September 2022
================================================================================================="""
import argparse

# import re

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # process command line arguments

    commandline = argparse.ArgumentParser(
        description='identify transcripts by diamond/blastx matches')
    commandline.add_argument('salmon',
                             help='salmon count file)',
                             type=argparse.FileType('r'))
    commandline.add_argument('blast',
                             help='tabular blast search result',
                             type=argparse.FileType('r'))
    # 'A:/mrg/Dropbox/colom/potato/trinity_out_dir.Trinity.fasta'
    # "/scratch/bell/mgribsko/220830uniref50/potato.trinity_uniref50_2.dmndblastx"

    cl = commandline.parse_args()

    print('\ntrinity_filterbysubject.py - filter by match\n')

    # read the index to count the total number of assembled sequences
    columns = cl.salmon.readline();
    counts = {}
    n_salmon = 0
    for line in cl.salmon:
        n_salmon += 1
        values = line.rstrip().split()
        count = 0
        for v in values[1:]:
            count += int(v)
        counts[values[0]] = count

    # look through the blastx result
    field = {'qid':   0, 'qlen': 1, 'qstart': 2, 'qstop': 3,
             'sid':   4, 'slen': 5, 'sstart': 6, 'sstop': 7,
             'allen': 8, 'pid': 9, 'score': 10, 'evalue': 11, 'stitle': 12}
    search_term = 'virus'

    n_component = 0
    n_match = 0
    n_result = 0
    component = {}
    for line in cl.blast:
        n_result += 1
        line = line.rstrip()
        # print(line)
        token = line.split('\t', maxsplit=12)
        title = token[field['stitle']]
        evalue = token[field['evalue']]
        length = int(token[field['qlen']])

        # construct the component name
        qtoken = token[field['qid']].split('_')
        qbaseid = '_'.join([qtoken[1], qtoken[2]])

        if qbaseid in component:
            component[qbaseid]['count'] += 1
        else:
            component[qbaseid] = {'count': 1, 'evalue': evalue, 'line': line, 'length': length}
            n_component += 1

        if evalue < component[qbaseid]['evalue']:
            component[qbaseid]['evalue'] = evalue
            component[qbaseid]['line'] = line
            component[qbaseid]['length'] = length

    # filter the best result for each component with the search term
    # for c in component:
    # for c in sorted(component, key=lambda c: component[c]['length']):

    n_nocount = 0
    for c in component:
        if not c in counts:
            print(f'no counts for component {c}')
            counts[c] = 0
            n_nocount += 1
    print(f'\nno counts: {n_nocount} set to zero')

    out = open('match.out','w')

    for c in sorted(component, key=lambda c: counts[c], reverse=True):
        # if subject title does not contain search term, skip
        if component[c]['line'].lower().find(search_term) == -1:
            continue

        n_match += 1
        out.write(f'{c}\t{counts[c]}\t{component[c]["line"]}\n')

    print(f'\nblast result: {n_result}\nsalmon: {n_salmon}\ncomponent: {n_component}\nmatch: {n_match}')

exit(0)

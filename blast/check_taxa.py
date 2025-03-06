"""#################################################################################################
After hand editing, check taxon lists to make sure all have the scientific name (Tax) as well as
NCBI numeric name (TaxID)

TODO add counters for good and bad files

26 February 2025        Michael Gribskov
#################################################################################################"""
from collections import defaultdict


def read_tsv(fname, ):
    """---------------------------------------------------------------------------------------------
    read the hand edited tsv file of taxa and taxid, tab separated

    field 0: count of matches in the blast search
    field 1: taxonomic name - continues to first all numerical field (may have mutliple terms)
    field 2: taxonomic id - the first all digit field after the


    :param fname: string    path to tsv data file
    :return: list of dict   data read from file
    ---------------------------------------------------------------------------------------------"""
    f = open(fname, 'r')
    out = defaultdict(lambda: {'species': set(), 'taxid': set()})
    for line in f:
        print(f'{line}')
        field = line.rstrip().split('\t')
        tax = field[1]
        this = out[tax]

        for i in range(2, len(field)):
            if field[i].isdigit():
                break
            this['species'].add(field[i])

        try:
            if i < len(field):
                this['taxid'].add(int(field[i]))
        except:
            print(f'missing taxid')

    f.close()
    return out


def blast_query_block(blast):
    """---------------------------------------------------------------------------------------------

    :param blast:       blast search result
    :return:
    ---------------------------------------------------------------------------------------------"""
    query_old = ''
    qinfo = {'query': '', 'genera': []}

    for line in blast:
        field = line.split('\t')
        query = field[0]
        taxbegin = field[12].find("Tax") + 4
        taxend = field[12].find("TaxID") - 1
        taxterm = field[12][taxbegin:taxend].split()

        if query == query_old:
            # add to genera
            qinfo['query'] = query
            qinfo['genera'].append(taxterm[0])
        else:
            yield qinfo
            query_old = qinfo['query']
            qinfo['query'] = query
            qinfo['genera'] = [taxterm[0]]

    return False


if __name__ == '__main__':

    blastfile = 'data/c16c31.trinity_uniref_1e-10.dmndblastx'
    goodfile = 'data/tax_good.txt'
    badfile = 'data/tax_bad.txt'

    good = read_tsv(goodfile)
    # manually add
    good['Potato'] = {'species': ['species'], 'taxid': [4081]}

    goodset = open('goodset.out', 'w')
    badset = open('badset.out', 'w')

    blast = open(blastfile, 'r')
    for qinfo in blast_query_block(blast):
        # check to see if any genus in qinfo['genera'] is in good list
        ok = False
        for genus in qinfo['genera']:
            if genus in good:
                ok = True
                break

        if ok:
            # print(f'GOOD {qinfo}')
            goodset.write(f"{qinfo['query']}\t{qinfo['genera']}\n")
        else:
            badset.write(f"{qinfo['query']}\t{qinfo['genera']}\n")
            # if 'root' in qinfo['genera'] or 'Eukaryota' in qinfo['genera']:
            #     print(f'BAD {qinfo}')

    goodset.close()
    badset.close()

    exit(0)

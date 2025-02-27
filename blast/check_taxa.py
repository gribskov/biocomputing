"""#################################################################################################
After hand editing, check taxon lists to make sure all have the scientific name (Tax) as well as
NCBI numeric name (TaxID)

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
    out = defaultdict(lambda: {'species':set(), 'taxid':set()})
    for line in f:
        print(f'{line}')
        field = line.rstrip().split('\t')
        tax = field[1]
        this = out[tax]

        for i in range(2,len(field)):
            if field[i].isdigit():
                break
            this['species'].add(field[i])

        this['taxid'].add(int(field[i]))


    return out

if __name__ == '__main__':

    blastfile = 'data/c16c31.trinity_uniref_1e-10.dmndblastx'
    goodfile = 'data/tax_good.txt'
    badfile = 'data/tax_bad.txt'

    good = read_tsv(goodfile)

    exit(0)

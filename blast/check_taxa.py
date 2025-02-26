"""#################################################################################################
After hand editing, check taxon lists to make sure all have the scientific name (Tax) as well as
NCBI numeric name (TaxID)

26 February 2025        Michael Gribskov
#################################################################################################"""

def read_tsv(fname):
    """---------------------------------------------------------------------------------------------
    read the hand edited tsv file of taxa and taxid

    field 0: count of matches in the blast search
    field 1: taxonomic name - continues to first all numerical field

    :param fname: string    path to tsv data file
    :return: list of dict   data read from file
    ---------------------------------------------------------------------------------------------"""
    f = open(fname, 'r')
    out = []
    for line in f:
        print(f'{line}')
        field = line.rstrip().split('\t')
        tax = field[1]
        for i in range(2,len(field)):
            if field[i].isdigit():
                break

            tax += f' {field[i]}'

        taxid = ' '.join(field[i:])
        if not taxid.isdigit():
            print(f'|{tax}|{taxid}|')


        out.append(line)

    return out

if __name__ == '__main__':

    blastfile = 'data/c16c31.trinity_uniref_1e-10.dmndblastx'
    goodfile = 'data/tax_good.txt'
    badfile = 'data/tax_bad.txt'

    good = read_tsv(goodfile)

    exit(0)

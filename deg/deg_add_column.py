"""=================================================================================================
add the potato functional annotation to the deg stats and counts (TSV)

DM_1-3_516_R44_potato.v6.1.working_models.func_anno.txt

Soltu.DM.01G001080.1	ARM repeat superfamily protein
Soltu.DM.01G001090.1	related to AP2
Soltu.DM.01G001100.1	Plant protein 1589 of unknown function (A_thal_3526) domain containing protein
Soltu.DM.01G001100.2	Plant protein 1589 of unknown function (A_thal_3526) domain containing protein
Soltu.DM.01G001100.3	Plant protein 1589 of unknown function (A_thal_3526) domain containing protein

Soltu.DM.03G026180.1 (IDs in TSV)

Michael Gribskov     05 December 2022
================================================================================================="""
import sys


def read_column(file, key_n, value_n, maxsplit):
    """---------------------------------------------------------------------------------------------
    read file and create a dict by splitting each line split_n times and using column key_n as the
    key and value_n as a value

    :param file: filehandle     open file for reading, must be splitable
    :param key_n: int           column for keys
    :param value_n: int         column for values
    :param maxsplit: int        number of times to split each line
    :return: dict
    ---------------------------------------------------------------------------------------------"""
    data = {}
    for line in file:
        field = line.rstrip().split('\t', maxsplit=maxsplit)
        data[field[key_n]] = field[value_n]

    return data


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print(f'deg_add_columns.py')

    annofile = sys.argv[1]
    anno = open(annofile, 'r')
    print(f'Annotation column in {annofile}')

    degfile = sys.argv[2]
    deg = open(degfile, 'r')
    print(f'DEG data in {degfile}')

    outfile = sys.argv[3]
    out = open(outfile, 'w')
    print(f'Output written to {outfile}')

    coldata = read_column(anno, 0, 1, 1)
    print(f'\n{len(coldata)} annotations read from {annofile}')
    anno.close()

    # add a column label to the first line
    gene_n = 0
    skipline = deg.readline().rstrip()
    skipline += '\tFunctional_Annotation\n'
    out.write(skipline)

    for line in deg:
        gene_n += 1
        try:
            (gene, data) = line.rstrip().split('\t', maxsplit=1)
        except ValueError:
            sys.stderr.write(f'too few fields to split in deg file:{line.rstrip()}\n')

        try:
            data += f'\t{coldata[gene]}'
        except IndexError:
            sys.stderr.write(f'gene {gene} not found in annotation\n')

        out.write(f'{gene}{data}\n')

    sys.stdout.write(f'{gene_n} genes read from {degfile} and written to {outfile}\n')
    out.close()
    deg.close()

    exit(0)

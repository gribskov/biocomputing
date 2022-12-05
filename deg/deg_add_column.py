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
    annofile = sys.argv[1]
    anno = open(annofile, 'r')
    print(f'Annotation column in {annofile}')

    coldata = read_column( anno, 0, 1, 1)

    exit(0)

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
import os


def read_subsplit(file, key_n, value_n, split='\t', subsplit=' '):
    """---------------------------------------------------------------------------------------------
    read file and create a dict by splitting each line on split. If field[key_n] is empty, no entry
    is created. Otherwise, the field is split on subsplit, and each token added to the dictionary.

    :param file: filehandle     open file for reading, must be splitable
    :param key_n: int           column for keys
    :param value_n: int         column for values
    :param split: string        first key to split on
    :param subsplit: string     second key, for splitting the desired field
    :return: dict
    ---------------------------------------------------------------------------------------------"""
    data = {}
    for line in file:
        try:
            field = line.rstrip('\n ').split(split)
            # field = line.rstrip().split(split)
        except IndexError:
            print('Cannot split {line')
            continue

        if len(field) <= value_n:
            continue

        if field[value_n]:
            indiv = field[value_n].split(subsplit)
            if field[key_n] in data:
                # existing gene
                data[field[key_n]] += indiv
            else:
                # new gene
                data[field[key_n]] = indiv

    return data


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


def construct_output_name(infile, inname, col_n, remove=['.tsv', '.txt'], header='auto',
                          dropfirst=False):
    """---------------------------------------------------------------------------------------------
    the output name is the input name minus any strings in del, plus the column name. If header is
    True, the column name in col_n position of the first line is used. otherwise the name is
    .anno_col_n. the header line is returned for re-use. If header is false, the first line is not
    read and hline is empty.

    header='auto'   read first line, parse column header from field col_n
    header='skip'   read first line, assign column as ,anno_{col_n}
    header='other'  read first line, assign column as ,anno_other
    header = ''     do not read first line, assign column as ,anno_{col_n}

    TODO no way to not read first line and give a name, hmm

    :param infile: fp          open file for reading
    :param inname: string      name of the input file
    :param col_n: int          column for data
    :param remove: list        strings to delete from input file name
    :param header: string      if 'false', header will be used for column title  (not 'False')
    :param dropfirst: logical  drop the first token on the header line (comment symbol or tag)
    :return outname, hline:    string, string name of output file, header line
    ---------------------------------------------------------------------------------------------"""
    # construct base of output file from input file
    # get base filename and open output file
    base = os.path.basename(inname)
    for suffix in remove:
        base = base.replace(suffix, '')

    hline = ''
    col = ''
    if header:
        # header is a non empty string, read the header line
        hline = infile.readline().rstrip()
        if header == 'auto' or header == '1':
            # split header line to get column label
            field = hline.split()
            if dropfirst:
                field = field[1:]
            col = field[col_n]
        elif header == 'skip' or header == '0':
            # do not try to find a column name, proceed to automatic name
            header = ''
        else:
            # header is a string to use as a column title
            col = header

    if not header:
        # header is empty, assign column using header (empty string is false)
        col = f'.anno_{col_n}'

    outname = f'{base}_{col}'

    return outname, hline


def write_append_deg(deg, out, coldata):
    """---------------------------------------------------------------------------------------------
    append to the column of interest to the end of the appropriate line of the deg file, separated
    by a tab.

    :param deg: fp          deg file (tsv) for reading
    :param out: fp          output file for writing
    :param coldata: dict    genes, labels
    :return: int            rec ords written
    ---------------------------------------------------------------------------------------------"""
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
        except KeyError:
            data += f'\t'
            sys.stderr.write(f'gene {gene} not found in annotation\n')

        out.write(f'{gene}{data}\n')

    return gene_n


def write_filtered_columm(deg, out, coldata):
    """---------------------------------------------------------------------------------------------
    write the data, removing any IDs that don't occur in deg

    :param deg: fp          deg data open for reading
    :param out: fp          output file open for writing
    :param coldata: dict    gene,label
    :return: int            lines written
    ---------------------------------------------------------------------------------------------"""
    # first make a list of the desired genes (first column of DEG file)
    accept = []
    for line in deg:
        field = line.split()
        if field[0] not in accept:
            accept.append(field[0])

    row_n = 0
    for gene in coldata:
        if gene in accept:
            for label in coldata[gene]:
                out.write(f'{gene}\t{label}\n')
                row_n += 1

    return row_n


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

    col_n = 7
    (outfile, header) = construct_output_name(anno, degfile, col_n, header='auto')
    out = open(outfile, 'w')
    print(f'Output written to {outfile}')

    # coldata = read_column(anno, 0, 1, 1)
    coldata = read_subsplit(anno, 2, col_n, split='\t', subsplit=' ')
    print(f'\n{len(coldata)} annotations read from {annofile}')
    anno.close()

    # gene_n = write_append_deg(deg, out, coldata)
    # sys.stdout.write(f'{gene_n} genes read from {degfile} and written to {outfile}\n')
    row_n = write_filtered_columm(deg, out, coldata)
    sys.stdout.write(f'{row_n} annotation labels read from {degfile} and written to {outfile}\n')

    out.close()
    deg.close()

    exit(0)

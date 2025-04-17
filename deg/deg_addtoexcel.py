"""=================================================================================================
Add information to DEG excel file (derived from DESeq2 output)
Blast search results

Michael Gribskov     15 April 2025
================================================================================================="""
# from openpyxl import load_workbook
import argparse
import textwrap as _textwrap
import pandas as pd


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """---------------------------------------------------------------------------------------------
    Custom formatter for argparse.  Less ugly breaking of lines.
    ---------------------------------------------------------------------------------------------"""

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        self._max_help_position = 30
        return _textwrap.wrap(text, 90)


def getopts():
    """---------------------------------------------------------------------------------------------
    Get command line options
    :return: namespace with options
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(description='Add annotation to gene list',
                                 formatter_class=CustomFormatter)
    cl.add_argument('--atype', type=str, default='excel',
                    help='annotation file type (excel)')
    cl.add_argument('--dtype', type=str, default='blast',
                    help='data source file type (blast)')
    cl.add_argument('annotation', type=str, help='input annotation file to augment')
    cl.add_argument('data', type=str, help='data file with annotation information')
    cl.add_argument('output', type=str, help='name of augmented output file')

    return cl.parse_args()


def blast_read(data, rownames, maxeval=1e-5, nhits=3):
    """---------------------------------------------------------------------------------------------
    read in the blast data to add to the annotation
    :param data: filehandle     open for reading
    :param rownames: list       names corresponding to query id in blast result
    :param maxeval: float       maximum evalue
    :param nhits: int           maximum number of hits to report for each query
    :return: dataframe
    ---------------------------------------------------------------------------------------------"""
    selected = {}
    for id, hitlist in blast_get_block(data):
        if id not in rownames:
            continue

        # split and assemble the desired data
        unsorted = []
        for hit in hitlist:
            field = hit.rstrip().split(maxsplit=12)
            # print(f'\t{field[12]}')

            evalue = float(field[11])
            hitstr = f'{evalue}  '
            hitstr += f'q:{field[2]}:{field[3]}/{field[1]}  '
            hitstr += f's:{field[6]}:{field[7]}/{field[5]}   '
            hitstr += f'{field[12][:field[12].find(' TaxID')].replace('Tax=', '')}'

            unsorted.append([evalue, hitstr])

        anno = ''
        annolist = []
        e_lowest = None
        n = 0
        for info in sorted(unsorted, key=lambda h: h[0]):
            if not e_lowest:
                e_lowest = info[0]
            n += 1
            if n > nhits:
                break
            anno += f'{info[1]}\n'
            annolist.append(f'{info[1]}\n')
        anno = anno[:-1]
        # print(f'{id}\n{anno}')
        selected[id] = [e_lowest, anno]

    df = pd.DataFrame.from_dict(selected, orient='index',columns=['evalue', 'blast_hits'])

    return df


def blast_get_block(data, level=3):
    """---------------------------------------------------------------------------------------------
    generator that returns a list of blast hit results that have a common trinity transcript id
    level 3 is the _g level

    :param data:
    :return:
    ---------------------------------------------------------------------------------------------"""
    hits = []
    jlevel = level + 1

    # get the first hit
    line = data.readline()
    full_id = line.split(maxsplit=1)[0]
    idfield = full_id.split('_')
    id = '_'.join(idfield[1:jlevel])
    hits.append(line.rstrip())
    id_old = id

    for line in data:
        full_id = line.split(maxsplit=1)[0]
        idfield = full_id.split('_')
        id = '_'.join(idfield[1:jlevel])

        if id == id_old:
            # same query as last line,  add to hits
            hits.append(line.rstrip())
        else:
            yield id, hits
            hits = [line.rstrip()]
            id_old = id

    if hits:
        yield id, hits


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    opts = getopts()

    data = open(opts.data, 'r')
    out = open(opts.output, 'w')

    # to get a list of worksheets in workbook
    # sheets = pd.ExcelFile(opts.annotation).sheet_names

    annodf = None
    if opts.atype == 'excel':
        # read the first sheet into a dataframe with row and column labels
        annodf = pd.read_excel(opts.annotation, sheet_name=0, index_col=0)
        print(f'Annotations read from {opts.annotation}: {annodf.shape}')
    #
    # colnames = list(annodf.columns)
    rownames = list(annodf.index)
    # print(f'columns: {colnames[0:4]}')
    # print(f'rows: {rownames[0:4]}')


    if opts.dtype == 'blast':
        blastdf = blast_read(data, rownames, maxeval=1e-5, nhits=3)
        print(f'Blast results read from {opts.data}: {blastdf.shape}')

    # print(blastdf.head())
    # merge blast results with annotatio
    annodf.merge(blastdf, how='left', left_index=True, right_index=True )
    annodf = annodf.merge(blastdf, how='outer', left_index=True, right_index=True)

    # print(annodf.head())
    # annodf.to_excel('test.xlsx')
    # print(annodf.iloc[1:4,67:])

    writer = pd.ExcelWriter('pandas_test.xlsx', engine='xlsxwriter')
    annodf.to_excel(writer, sheet_name='Annotation', index=True)
    workbook = writer.book
    worksheet = writer.sheets['Annotation']

    # Add a text wrap format and format the column with the multi-line blast info
    text_wrap_format = workbook.add_format({'text_wrap': True})
    worksheet.set_column(69, 69, 200, text_wrap_format)

    writer.close()
    data.close()
    out.close()
    exit(0)

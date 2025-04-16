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
    for id,hitlist in blast_get_block(data):
        if id not in rownames:
            continue

        print(f'{id}')
        for hit in hitlist:
            print(f'\t{hit}')

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

    colnames = list(annodf.columns)
    rownames = list(annodf.index)
    print(f'columns: {colnames[0:4]}')
    print(f'rows: {rownames[0:4]}')

    if opts.dtype == 'blast':
        blastdf = blast_read(data, rownames, maxeval=1e-5, nhits=3)

    data.close()
    out.close()
    exit(0)

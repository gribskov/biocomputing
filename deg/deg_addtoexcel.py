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


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    opts = getopts()
    # xlsx = 'A:/mrg/Dropbox/colom/potato/2411_C16C31/C16C31.lfc.xlsx'
    # wb = load_workbook(xlsx)
    #
    # # grab the active worksheet
    # ws = wb.active
    # # sheet_ranges = wb['range names']
    # # print(sheet_ranges['D18'].value)
    #
    # # Data can be assigned directly to cells
    # # ws['A1'] = 42
    #
    # for i in range(10):
    #     value = ws(f'A{i}')
    #     print(f'{i}: {value}')

    # anno = open(opts.annotation, 'r')
    data = open(opts.data, 'r')
    out = open(opts.output, 'w')

    # to get a lest of worksheets in workbook
    # sheets = pd.ExcelFile(opts.annotation).sheet_names
    if opts.atype == 'excel':
        # read the first sheet into a dataframe with row and column labels
        annodf = pd.read_excel(opts.annotation, sheet_name=0, index_col=0)
    # sheet = workbook[sheets[0]]
    # print(sheet, sheets)

    print(annodf.iloc[0:6, 0:3])
    colnames = annodf.columns.tolist()
    colnames = list(annodf.columns)
    rownames = list(annodf.index)
    print(f'columns: {colnames[0:4]}')
    print(f'rows: {rownames[0:4]}')

    # anno.close()
    data.close()
    out.close()
    exit(0)

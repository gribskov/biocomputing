"""=================================================================================================
Add information to DEG excel file (derived from DESeq2 output)
Blast search results

Michael Gribskov     15 April 2025
================================================================================================="""
# from openpyxl import load_workbook
import pandas as pd

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    xlsx = 'A:/mrg/Dropbox/colom/potato/2411_C16C31/C16C31.lfc.xlsx'
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

    sheets = pd.ExcelFile(xlsx).sheet_names
    df = pd.read_excel(xlsx, sheet_name=sheets,)

    print(sheets)
    print(df.head())
    print(df.columns.tolist())
    print(df.iloc[0:6,0:3])


    exit(0)

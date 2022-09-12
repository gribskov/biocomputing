"""=================================================================================================


Michael Gribskov     11 April 2021
================================================================================================="""
import codewars_test as test

t = {'RR':'R', 'RB':'G', 'RG':'B',
     'BR':'G', 'BB':'B', 'BG':'R',
     'GR':'B', 'GB':'R', 'GG':'G'}


def triangle(row):
    while len(row) > 1:
        rowlist = [t[row[i:i + 2]] for i in range(len(row) - 1)]
        row = ''.join(rowlist)

    return row


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    test.assert_equals(triangle('GB'), 'R')
    test.assert_equals(triangle('RRR'), 'R')
    test.assert_equals(triangle('RGBG'), 'B')
    test.assert_equals(triangle('RBRGBRB'), 'G')
    test.assert_equals(triangle('RBRGBRBGGRRRBGBBBGG'), 'G')
    test.assert_equals(triangle('B'), 'B')

    exit(0)

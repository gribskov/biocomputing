"""=================================================================================================
Rule 30 is a one-dimensional binary cellular automaton. You can have some information here: https://en.wikipedia.org/wiki/Rule_30

Complete the function that takes as input an array of 0s and 1s and a non-negative integer n that represents the number of iterations. This function has to perform the nth iteration of Rule 30 with the given input.

The rule to derive a cell from itself and its neigbour is:

Current cell	000	001	010	011	100	101	110	111
New cell	0	1	1	1	1	0	0	0
As you can see the new state of a certain cell depends on his neighborhood. In Current cell you have the nth cell with its left and right neighbor, for example the first configuration is 000:

left neighbor = 0
current cell = 0
right neighbor = 0
The result for the current cell is 0, as reported in the New cell row.

You also have to pay attention to the following things:

the borders of the list are always 0
you have to return an array of 0 and 1
Here a small example step by step, starting from the list [1] and iterating 5 times:

Michael Gribskov     09 June 2025
================================================================================================="""

def rule30(cell, n):
    xform = {'000': '0', '001': '1', '010': '1', '011': '1', '100': '1', '101': '0', '110': '0', '111': '0'}
    cell = ''.join(str(c) for c in cell)
    i = 0
    while i < n:
        cell = '00' + cell + '00'
        # print(cell)
        l = len(cell)

        cellnew = ''
        for pos in range(1, l - 1):
            env = cell[pos - 1:pos + 2]
            cellnew += xform[env]

        cell = cellnew
        i += 1

        # print('\t', cellnew)

    return [int(c) for c in cell]
# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
print(rule30([1], 1), [1, 1, 1])
print(rule30([1], 2), [1, 1, 0, 0, 1])

exit(0)

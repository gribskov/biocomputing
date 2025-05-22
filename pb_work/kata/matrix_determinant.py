"""=================================================================================================
Write a function that accepts a square matrix (N x N 2D array) and returns the determinant of the matrix.

How to take the determinant of a matrix -- it is simplest to start with the smallest cases:

A 1x1 matrix |a| has determinant a.

A 2x2 matrix [ [a, b], [c, d] ] or

|a  b|
|c  d|
has determinant: a*d - b*c.

The determinant of an n x n sized matrix is calculated by reducing the problem to the calculation of the determinants of n matrices ofn-1 x n-1 size.

For the 3x3 case, [ [a, b, c], [d, e, f], [g, h, i] ] or

|a b c|
|d e f|
|g h i|
the determinant is: a * det(a_minor) - b * det(b_minor) + c * det(c_minor) where det(a_minor) refers to taking the determinant of the 2x2 matrix created by crossing out the row and column in which the element a occurs:

|- - -|
|- e f|
|- h i|
Note the alternation of signs.

The determinant of larger matrices are calculated analogously, e.g. if M is a 4x4 matrix with first row [a, b, c, d], then:

det(M) = a * det(a_minor) - b * det(b_minor) + c * det(c_minor) - d * det(d_minor)

Michael Gribskov     20 May 2025
================================================================================================="""


def do_test(mat, expected):
    print(f'matrix determinant {mat}:', end=' ')
    val = determinant(mat)
    print(f'result={val} / expected={expected} ', end='\t')
    status = 'failed'
    if val == expected:
        status = 'passed'
    print(status)

    return


def determinant(mat):
    """---------------------------------------------------------------------------------------------
    calculate the determinant of a matrix recursively
    :param mat: list of list    input matrix
    :return: int                value of determinant
    ---------------------------------------------------------------------------------------------"""
    size = len(mat)
    if size == 1:
        return mat[0][0]
    # if size == 2:
    #     return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1]

    sign = 1
    det = 0
    for i in range(size):
        det += sign * mat[0][i] * determinant(minor(i, mat))
        sign *= -1

    return det


def minor(col, mat):
    """---------------------------------------------------------------------------------------------
    return the minor matrix for the specified column
    :param pos: int             col for minor generation
    :param mat: list of list    matrix
    :return: list of list       minor matrix
    ---------------------------------------------------------------------------------------------"""
    size = len(mat)
    minor = [[] for _ in range(size - 1)]
    r = 0
    for row in mat[1:]:
        for c in range(size):
            if c == col:
                continue
            minor[r].append(row[c])
        r += 1

    return minor


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
m1 = [[4, 6], [3, 8]]
m5 = [[2, 4, 2], [3, 1, 1], [1, 2, 0]]

do_test([[5]], 5)
do_test(m1, 14)
do_test(m5, 10)

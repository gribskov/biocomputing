"""=================================================================================================
Solution:
Generate all possible permutations of integers 1-n
evaluate each for "visibility" from forward and backward direction
based on the clues, make a list of the allowed permutations in each row and column
iterate over:
    generate all possile solutions from constraints
    update constrants based on current solution
until only one solution remains

Michael Gribskov     04 April 2021
================================================================================================="""
from itertools import permutations
import codewars_test as test


def apply_solution(size, solution, constraints):
    """
    Remove constraints that are inconsistent with the solution

    :param solution:
    :param constraints:
    :return:
    """
    n_new = 0
    new = {'row':[[] for _ in range(size)],
           'col':[[] for _ in range(size)]}

    # use possible  solutions to get improved row constraints
    for r in range(size):
        for possibility in constraints['row'][r]:
            possible = True
            for c in range(size):
                if possibility[c] not in solution[r][c]:
                    possible = False
                    break
            if possible:
                new['row'][r].append(possibility)
                n_new += 1

    # use possible  solutions to get improved col constraints
    for c in range(size):
        for possibility in constraints['col'][c]:
            possible = True
            for r in range(size):
                if possibility[r] not in solution[r][c]:
                    possible = False
                    break
            if possible:
                new['col'][c].append(possibility)
                n_new += 1

    return n_new, new


def apply_constraints(size, constraints):
    """
    The constraints give the permutations that are currently allowed for each column (top->down) and
    row (left->right).  the intersection of the of the row and column constraint in each cell is the
    set of possible solutions for the cell.

    :param constraints:
    :return: solution
    """
    solution = [[None for _ in range(size)] for _ in range(size)]
    for r in range(size):
        for c in range(size):
            solution[r][c] = intersect(constraints, r, c)

    return solution


def intersect(constraints, r, c):
    """
    Get the intersection of the row and column constraints for row r and column c
    :param constraints:
    :param r:
    :param c:
    :return: list, possible solutions
    """
    rowset = set()
    colset = set()
    arow = constraints['row'][r]
    for possibility in constraints['row'][r]:
        rowset.add(possibility[c])
    for possibility in constraints['col'][c]:
        colset.add(possibility[r])

    rowset.intersection_update(colset)
    return list(rowset)


def setup_constraints(clues, matrix):
    """
    based on the clues, find the allowed permutations for each tow and column
    :param clues:
    :param matrix:
    :return: int:n_constraints, dict: constraints
    """
    clues = list(clues)
    size = len(clues)
    quarter = size // 4

    constraints = {'row':[None for _ in range(quarter)],
                   'col':[None for _ in range(quarter)]}

    # get the constraints on the rows and columns from the clues
    # each constraint is two number reflecting the number of visible skyscrapers from the two
    # directions: top->down and left->right
    for n in range(quarter):
        constraints['col'][n] = [clues[n], clues[size - quarter - n - 1]]
        constraints['row'][n] = [clues[size - 1 - n], clues[n + quarter]]

    # replace constraints with the sets of solutions possible in the rows and columns
    n_constraints = 0
    for k in ('row', 'col'):
        for pos in range(len(constraints[k])):
            allowed = matrix[constraints[k][pos][0]][constraints[k][pos][1]]
            # constraints[k][pos] = matrix[constraints[k][pos][0]][constraints[k][pos][1]]
            constraints[k][pos] = allowed
            n_constraints += len(allowed)

    return n_constraints, constraints


def visible(order):
    """
    Count the number of scyscrapers visible for a permutation
    :param order:
    :return:
    """
    visible = 0
    last = 0
    for i in order:
        if i > last:
            visible += 1
            last = i

    return visible


def permute(list):
    """
    Generate a list of permutation of all integers in list.  For each permutaiton calculate how many
    skyscrapers can be seen in the forward and backward direction.  Store indexed by forward and
    backward.
    Zero indicates 'no clue' so all permutations are visible allowed.
    :param list:
    :return:
    """
    n = len(list)
    matrix = [[[] for _ in range(n + 1)] for _ in range(n + 1)]
    p = permutations(list)
    for x in p:
        forward = visible(x)
        backward = visible(x[::-1])
        matrix[forward][backward].append(x)
        matrix[0][backward].append(x)
        matrix[forward][0].append(x)
        matrix[0][0].append(x)

    return matrix


def solve_puzzle(clues):
    size = 4
    solution = []
    values = [x for x in range(1, size + 1)]
    allowed = permute(values)
    possible, constraints = setup_constraints(clues, allowed)
    while possible > size * 2:
        solution = apply_constraints(size, constraints)
        possible, constraints = apply_solution(size, solution, constraints)

    # final = tuple((solution[r][c][0] for c in range(size)) for r in range(size))
    final = []
    for r in range(size):
        for c in range(size):
            solution[r][c] = solution[r][c][0]
        t = tuple(solution[r])
        final.append(t)

return tuple(final)


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    clues = (
        (2, 2, 1, 3,
         2, 2, 3, 1,
         1, 2, 2, 3,
         3, 2, 1, 3),
        (0, 0, 1, 2,
         0, 2, 0, 0,
         0, 3, 0, 0,
         0, 1, 0, 0)
        )

    outcomes = (
        ((1, 3, 4, 2),
         (4, 2, 1, 3),
         (3, 4, 2, 1),
         (2, 1, 3, 4)),
        ((2, 1, 4, 3),
         (3, 4, 1, 2),
         (4, 2, 3, 1),
         (1, 3, 2, 4))
        )

    test.describe("4 by 4 skyscrapers")
    test.it("should pass all the tests provided")

    test.assert_equals(solve_puzzle(clues[0]), outcomes[0])
    test.assert_equals(solve_puzzle(clues[1]), outcomes[1])

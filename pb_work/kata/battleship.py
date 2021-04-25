import codewars_test as test

directions = [(-1, 0),(1, 0), (0, -1), (0, 1)]
diagonal = [(-1, -1), (-1, 1), (1, -1), (1, 1)]

def traceship(field, row, col):
    ship = []
    if row < 0 or row > 9 or col < 0 or col > 9:
        return ship

    if not field[row][col]:
        return ship

    diag = False
    for offset in diagonal:
        x = row+offset[0]
        y = col+offset[1]
        if x < 0 or x > 9 or y < 0 or y > 9:
            continue
        if field[x][y]:
            diag = True
            break
    if diag:
        return ship

    field[row][col] = 0
    ship = [(row, col)]
    for offset in directions:
        ship += traceship(field, row + offset[0], col + offset[1])

    return ship


def validate_battlefield(field):
    shiplist = []
    for row in range(len(field)):
        for col in range(len(field[row])):
            ship = traceship(field, row, col)
            if ship:
                shiplist.append(ship)

    count = [0 for _ in range(5)]
    for ship in shiplist:
        try:
            count[len(ship)] += 1
        except IndexError:
            return False

    if count[4] != 1 or count[3] != 2 or count[2] != 3 or count[1] != 4:
        return False

    print(count)
    print(shiplist)
    return True


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    battleField = [[1, 0, 0, 0, 0, 1, 1, 0, 0, 0],
                   [1, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                   [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
                   [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

    test.assert_equals(validate_battlefield(battleField), True, "Yep! Seems alright",
                       "Nope, something is wrong!");
    exit(0)

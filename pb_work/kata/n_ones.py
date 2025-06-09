"""=================================================================================================
Background
Most numbers contain the digit 1; some do not. In fact, a whopping 53% of integers below 10^7 (ten million) contain the digit 1! As we consider larger and larger scales, the percentage of integers that contain the digit 1 tends towards 100%.

Problem Description
In this kata, you are required to build a function that receives a positive integer argument, n, and returns the n-th positive integer containing the digit 1, the first of which is 1.

The following numbers contain the digit 1: 10, 21, 1024, 1111111.

The following numbers do not contain the digit 1: 42, 666, 2048, 234567890.

The first 10 integers containing the digit 1 are 1, 10, 11, 12, 13, 14, 15, 16, 17 and 18.

Input Constraints
Fixed test cases: 1 ≤ n ≤ 100 (one hundred)

Small test cases: 1 ≤ n ≤ 10^5 (one hundred thousand)

Medium test cases: 1 ≤ n ≤ 10^10 (ten billion)

Large test cases: 1 ≤ n ≤ 10^15 (one quadrillion)

Michael Gribskov     09 June 2025
================================================================================================="""

def cell_num_containing_ones(n):
    inc = 1
    next = 10
    count = 0
    while inc:
        count += inc
        if count == n:
            return count

        if count + inc * 10 > n:
            # can't go to bigger poser stay at current inc
            pass
        elif count + inc > n:
            inc //= 10


    return count


    digits = len(s_start) - 1
    mult = 10 ** digits
    n = 0
    for digit in s_start:
        if digit == '1':
            n += mult
        mult //= 10

    return n

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
print(cell_num_containing_ones(1), 1, 'n = 1')
print(cell_num_containing_ones(2), 10, 'n = 2')
print(cell_num_containing_ones(3), 11, 'n = 3')
print(cell_num_containing_ones(4), 12, 'n = 4')
print(cell_num_containing_ones(5), 13, 'n = 5')
print(cell_num_containing_ones(10), 18, 'n = 10')
print(cell_num_containing_ones(19), 91, 'n = 19')
print(cell_num_containing_ones(20), 100, 'n = 20')
print(cell_num_containing_ones(50), 130, 'n = 50')
print(cell_num_containing_ones(100), 180, 'n = 100')

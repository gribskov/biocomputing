"""=================================================================================================


Michael Gribskov     09 June 2025
================================================================================================="""
def sum_of_intervals(ranges):

    current = None
    result = 0
    for r in sorted(ranges, key=lambda r:r[0]):
        if current:
            if r[0] < current[1]:
                current[1] = max(current[1], r[1])
            else:
                result += current[1] - current[0]
                current = list(r)
            continue
        current = list(r)

    result += current[1] - current[0]

    return result

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
# print(sum_of_intervals([(1, 5)]), 4)
# print(sum_of_intervals([(1, 5), (6, 10)]), 8)
# print(sum_of_intervals([(1, 5), (1, 5)]), 4)
# print(sum_of_intervals([(1, 4), (7, 10), (3, 5)]), 7)
print(sum_of_intervals([(-347, 328), (-444, -128), (132, 186), (-372, 133)]), 772)

print(sum_of_intervals([(-1_000_000_000, 1_000_000_000)]), 2_000_000_000)
print(sum_of_intervals([(0, 20), (-100_000_000, 10), (30, 40)]), 100_000_030)

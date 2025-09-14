"""=================================================================================================


Michael Gribskov     11 June 2025
================================================================================================="""


def rr(n):
    circle = [x + 1 for x in range(n)]
    result = []

    for r in range(n - 1):
        schedule = []

        for i in range(n // 2):
            schedule.append((circle[i], circle[n - i - 1]))
        save = circle[n - 1]
        for i in range(n - 1, 1, -1):
            circle[i] = circle[i - 1]
        circle[1] = save

        result.append(schedule)

    return result


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------

print(rr(4))
print()
print(rr(6))

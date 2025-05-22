"""=================================================================================================


Michael Gribskov     22 May 2025
================================================================================================="""
import sys
from collections import defaultdict


def factor(n):
    """
    report prime factors of n
    :param n:
    :return: dict       key=factor, value = count
    """
    f = defaultdict(int)
    if n == 0: return f

    while not n % 2:
        f[2] += 1
        n //= 2

    for i in range(3, n + 1, 2):
        while not n % i:
            f[i] += 1
            n //= i
        i += 2
        if i > n:
            break

    return f


def count_factors(num, factors):
    """
    for each factor in the base , count how many occur in num
    :param num: int                 number to generate factorial
    :param factors: defaultdict     count of prime factors in base
    :return: defayultdict           count of base factors in num
    """
    count = defaultdict(int)
    for f in factors:
        base = f
        while base <= num:
            count[f] += num // base
            base *= f
    return count


def trailing_zeros(num, base):
    """
    a zero is generated for every set of prime factors found in the base. The
    number of zeros is limited by the total count / number in base
    :param num: int         number to generate factorial
    :param base: int        base for representing the factorial
    :return: int            number of trailing zeros
    """
    f = factor(base)
    c = count_factors(num, f)
    zeros = num
    for i in c:
        zeros = min(zeros, c[i] // f[i])

    return zeros


# main

test = [(15, 10, 3), (7, 21, 1), (15, 12, 5)]

for t in test:
    result = trailing_zeros(t[0], t[1])
    if result == t[2]:
        print(f'{t} passed')
    else:
        print(f'{t} failed {result}')

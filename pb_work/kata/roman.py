"""=================================================================================================


Michael Gribskov     11 June 2025
================================================================================================="""
from sys import prefix

r2a = { 'M': 1000,
        'CM': 900,
        'D': 500,
        'CD': 400,
        'C': 100,
        'XC': 90,
        'L': 50,
        'XL': 40,
        'X': 10,
        'IX': 9,
        'V': 5,
        'IV': 4,
        'I': 1, }

a2r = {r2a[r]:r for r in r2a}

def from_roman(roman):
    arabic = 0
    while roman:
        for prefix in sorted(r2a, key=lambda r:r2a[r], reverse=True):
            if roman.startswith(prefix):
                arabic += r2a[prefix]
                roman = roman[len(prefix):]
                break

    return arabic

def to_roman(arabic):
    roman = ''
    while arabic > 0:
        for rval in sorted(a2r, reverse=True):
            if arabic >= rval:
                roman += a2r[rval]
                arabic -= rval
                break

    return roman


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
print(from_roman('MMXXVI'))
print(to_roman(2026))

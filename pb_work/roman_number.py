def romanToint(str):
    # str = input("")
    nums = 0
    # for i in range(len(str)-1):
    i = 0
    while i < len(str) -  1:
        if roman(str[i]) < roman(str[i+1]):
            nums += roman(str[i+1]) - roman(str[i])
            i += 2
        else:
            nums += roman(str[i])
            i += 1

    # make sure to add the last digit, but not if it was added above
    if i < len(str):
        nums += roman(str[i])

    return nums

def roman(str):
    # dictionary = {'I':'1', 'V':'5', 'X':'10', 'L':'50', 'C':'100', 'D':'500', 'M':'1000'}
    # i = int(dictionary.values(str))
    dictionary = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    i = dictionary[str]
    return i

def romanII(str):
    """
    Convert roman numeral to string
    :param str: string with latin alphabet nuber (roman number)
    :return: int, arabic number
    """
    special = {'IV':4, 'IX':9, 'XL':40, 'XC':90, 'CD':400, 'CM':900}  # the only valid digraphs
    standard = {'I':1, 'V':5, 'X':10, 'L':50, 'C':100, 'D':500, 'M':1000 }

    number = 0

    # first convert the digraphs, digraphs can occur only once
    for di in special:
        if di in str:
            number += special[di]
            str = str.replace(di, '')

    for std in standard:
        number += str.count(std) * standard[std]

    return number



if __name__ == '__main__':
    test = ['III', 'IV', 'XXIX', 'LXVI', 'MCMLXVIII']
    for rom in test:
        print(f'roman {rom}\t{romanII(rom)}')

exit(0)


"""#################################################################################################
add_big_numbers.py


2025-05-23 gribskov
#################################################################################################"""
import sys


def tae(a, b, result):
    answer = add(a, b)
    print(f'add({a}, {b}) = {answer}  (result)', end=' ')
    status = 'fail'
    if answer == result:
        status = 'pass'
    print(f'{status}')
    return


def setup():
    """
    set up addition table
    :return:
    """
    plus = []
    for row in range(10):
        plus.append([])
        for col in range(10):
            s = row + col
            carry = 0
            if s > 9:
                carry = 1
                s = s - 10
            plus[row].append((s, carry))
    return plus


def add(a, b):
    # plus = setup()

    # ra = a.reverse()
    # rb = b.reverse()
    # make sure a is longer than b (or equal)
    if len(a) < len(b): a, b = b, a
    la = [int(c) for c in a[::-1]]
    lb = [int(c) for c in b[::-1]]

    result = ''
    carry = 0
    pa = len(la)
    pb = len(lb)
    for p in range(pa):
        if p >= pb:
            bval = 0
        else:
            bval = lb[p]
        # print(f'a[{pa}],b[{pb}]')
        s = la[p] + bval + carry
        carry = 0
        if s > 9:
            carry = 1
            s -= 10
        result += str(s)

    if carry: result += '1'



    return result[::-1];


####################################################################################################
# Main
####################################################################################################
if __name__ == '__main__':
    tae("1", "1", "2");
    tae("123", "456", "579");
    tae("888", "222", "1110");
    tae("1372", "69", "1441");
    tae("12", "456", "468");
    tae("101", "100", "201");
    tae('63829983432984289347293874', '90938498237058927340892374089', "91002328220491911630239667963")

exit(0)

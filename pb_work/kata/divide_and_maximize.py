"""=================================================================================================
Your friend has a list of k numbers: [a1, a2, a3, ... ak].

He is allowed to do an operation which consists of three steps:

select two numbers: ai and aj (ai % 2 = 0)
replace ai with ai / 2
replace aj with aj * 2
Help him to find the maximum sum of list elements that is possible to achieve by using this operation (possibly multiple times).
Return this sum modulo 1_000_000_007, because it can be quite big.

Input
List of k elements: [a1, a2, a3, ..., ak]; k < 10**4
All numbers are positive and smaller than 10**9

Output
Maximum sum after some operations (modulo 1_000_000_007)

Michael Gribskov     19 May 2025
================================================================================================="""


def do_test(ns, expected):
    copy = ns[:]
    val = divide_and_multiply(ns)
    print(f'divide_and_multiply {copy}: result={val} / expected={expected} ', end='\t')
    status = 'failed'
    if val == expected:
        status = 'passed'
    print(status)
    return


def parse(raw):
    """---------------------------------------------------------------------------------------------

    :param raw:
    :return:
    ---------------------------------------------------------------------------------------------"""
    div = []
    mult = []
    for i in range(len(raw)):
        if (not raw[i] % 2):
            div.append(i)
        mult.append(i)

    return div, mult


def divide_and_multiply_stack(ns):
    if len(ns) < 2:
        # must have two or more digits
        return sum(ns)

    # div = None
    # mult = None
    stack = []
    seen = []
    ns.sort()
    stack.append(ns)
    maxsum = 0
    while stack:
        current = stack.pop()
        maxsum = max(maxsum, sum(current))
        print(f'stack:{len(stack)}\t{current}\t{maxsum}')

        if current in seen:
            continue
        else:
            seen.append(current)
            div, mult = parse(current)
            if len(div) == 0:
                # cannot apply divide and multiply operation
                continue

            # add mult div combos to stack
            for d in div:
                tmp = current[:]
                tmp[d] = current[d] // 2
                save = tmp[:]
                for m in mult:
                    if m == d:
                        continue
                    tmp[m] = current[m] * 2
                    tmp.sort()
                    if tmp.sort() not in seen:
                        stack.append(tmp)
                    tmp = save

    return maxsum


def divide_and_multiply(ns):
    """---------------------------------------------------------------------------------------------
    the largest sum is the largest non-2 factor times the product of all the 2-factors, plus the
    other non-2 factors

    :param ns: list of int      starting list of numbers to maximize
    :return: int                largest sum of positions in ns mod 1_000_000_007
    ---------------------------------------------------------------------------------------------"""
    multiplier = 1
    for i in range(len(ns)):
        while not ns[i] % 2:
            ns[i] = ns[i] // 2
            multiplier *= 2

    ns.sort()
    ns[-1] *= multiplier
    return sum(ns) % 1_000_000_007


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
print("divide_and_multiply")

do_test([10], 10)
do_test([3, 4], 13)
do_test([6, 4, 2], 50)
do_test([6, 4, 5], 44)
do_test([1, 2, 3, 4, 5], 46)
do_test([8] * 15, 371842558)
do_test([14, 9, 26], 68)

exit(0)

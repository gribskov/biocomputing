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
stack = []
seen = []


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
        if (raw[i] % 2):
            div.append(i)
        mult.append(i)

    def divide_and_multiply(ns):

    # Let the force be with you, warrior
    stack.append(ns)
    total = 0
    while stack:
        current = stack.pop()
        current.sort()
        test = sum(current)
        if test > total:
            total = test

        if current not in seen:
            seen.append(current)
            even, odd = parse(current)

            #TODO add mult div combos to stack


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

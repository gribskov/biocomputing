"""=================================================================================================


    
================================================================================================="""


def factorial(n):
    """---------------------------------------------------------------------------------------------
    recursive factorial function
    factorial(n) = n * factorial(n-1)
    ---------------------------------------------------------------------------------------------"""
    if n > 1:
        f = n * factorial(n - 1)
        return f

    return 1


def factorialStack(n):
    """----------------------------------------------------------------------------------
    factorial using stack
    ----------------------------------------------------------------------------------"""
    # create stack
    stack = []
    while n > 0:
        stack.append(n)
        n -= 1

    # multiply all stack values
    f = 1
    while stack:
        f *= stack.pop()

    return f

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
import timeit

for i in range(50):

    setup = 'from __main__ import factorial, factorialStack'

    recurse = 'factorial({})'.format(i)
    print('i={}\trecursive: {:.5f}'.format(
        i, timeit.timeit(stmt=recurse, setup=setup,number=100000)), end='\t')

    stack = 'factorialStack({})'.format(i)
    print('stack: {:.5f}'.format(timeit.timeit(stmt=stack, setup=setup, number=100000)))
"""=================================================================================================
Calculates and prints the values of the fibonacci series, F, up to max_n terms
F(n) = F(n-1) + f(n-2)

Michael Gribskov     18 January 2021
================================================================================================="""

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # initialize n_max, the maximum term to calculate
    n_max = 6

    # initialize the first two values, f and f_old
    f_old = 0
    f = 1

    # print the first two terms
    print('F(', 0, ') =', f_old)
    print('F(', 1, ') =', f)

    # for each successive term, up to n_max, print the new term
    n = 2
    while n < n_max:
        f_new = f + f_old
        print('F(', n, ') =', f_new)

        f_old = f
        f = f_new
        n += 1

    exit(0)

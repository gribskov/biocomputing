"""=================================================================================================


Michael Gribskov     26 February 2021
================================================================================================="""

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    star = '*'
    maxlen = 5

    nstar = int(maxlen * (maxlen + 3) / 2)
    star_list = [star] * nstar
    pos = -1
    for row in range(maxlen, 0, -1):
        pos += row + 1
        star_list[pos] = '\n'

    print(''.join(star_list), end='')
    print(':')

    for row in range(maxlen, 0, -1):
        print(star * row)
    print(':')

    stars = star * maxlen
    for row in range(maxlen):
        print(stars[row:maxlen])
    print(':')

    print('\n'.join([ row * star for row in range(maxlen,0,-1)]))
    print(':')

    # def nstr(n):
    #     return star * n
    #
    # result = map(lambda r:star*r,range(maxlen,0,-1))
    # print('map')
    print('\n'.join(list(map(lambda r:star*r,range(maxlen,0,-1)))))
    print(':')
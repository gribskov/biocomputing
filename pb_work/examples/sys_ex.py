"""=================================================================================================
system examples

Michael Gribskov     24 March 2021
================================================================================================="""
import sys

# command line
n = 0
for arg in sys.argv:
    print('argument {} => {}'.format(n, sys.argv[n]))
    n += 1


exit(0)

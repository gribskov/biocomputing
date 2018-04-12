"""=================================================================================================
regular expression examples

12 April 2018
================================================================================================="""
import re

s = 'To be or not to be, that is the question\n'
if re.search('A', s):
    print('A matches')
if re.search('t', s):
    print('t matches')
if re.search('n$', s):
    print('n found at end of line')
if re.search('n\n', s):
    print('n newline matches')

exit(0)
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
    print('n newline matches\n')

if re.match('to', s):
    print('to matches with match')
if re.search('to', s):
    print('to matches with search')
if re.match('To', s):
    print('To matches with match')
if re.match('[Tt]o', s):
    print('[Tt]o matches with match')
if re.search('^to', s):
    print('^to matches with search')
if re.search('^[Tt]o', s):
    print('^[Tt]o matches with search')
if re.search('^[Tt]o', s):
    print('^[Tt]o matches with search\n')

s = 'ACGAGTCTCGTAGTCGTAAAAATGTGTGT'

if re.search('N*', s):
    print('N* matches')
if re.search('N+', s):
    print('N+ matches')
if re.search('CGTA?', s):
    print('CGT followed by 0/1 A matches')
if re.search('A{2,5}', s):
    print('2 to 5 A found\n')

s = '12345\tACGAGTCTCGTAGTCGTAAAAATGTGTGT'
print('string:', s)
if re.search('\d+ \w+', s):
    print('\d+ \w+ matches')
if re.search('\d+\t+\w+', s):
    print('\d+\\t+\w+ matches')
if re.search('\d+\s+\w+', s):
    print('\d+\s+\w+ matches\n')

s1 = '12345\tACGAGTCTCGTAGTCGTAAAAATGTGTGT'
print('string:', s1)
if re.search('\d+\W?[ACGT]+', s1):
    print('\d+\W?[ACGT]+ matches')

s2 = s1 + 'N'
print('string:', s2)
if re.search('\d+\W?[ACGT]+', s2):
    print('\d+\W?[ACGT]+ matches')
if re.search('\d+\W?[^N]+', s2):
    print('\d+\W?[^N]+ matches')
if re.search('\d+\W?[^N]+$', s2):
    print('\d+\W?[^N]+$ matches')
print()

s = 'To be or not to be, that is the question'
match = re.search('[Tt]o', s)
print('string={}'.format(s))
start = match.start()
end = match.end()
print('regex=[Tt]o    matches {} at {} to {}'.format(s[start:end + 1], start, end))
matches = re.findall('[Tt]o', s)
for match in matches:
    print(match)

s = 'To be or not to be, that is the question'
match = re.search('([Tt]o).+(to)', s)
print('groups', match.group(0), '|', match.group(1) )

exit(0)

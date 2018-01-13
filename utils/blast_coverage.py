"""=====================================================================================================================
blast_coverage

Determine which regions of a query sequence are covered by subjects in a database search.
Blast search should be run with -outfmt 7 (tabular output)

usage
blast_coverage.py <blast_search_file>
====================================================================================================================="""
import sys

# open the input file
blast = 0
try:
    blast = open(sys.argv[1], 'r')
except:
    print('Usage: blast_coverage.py <blast_search>')
    print('    unable to open input search file ({})'.format(sys.argv[1]))
    exit(1)

# skip the comment block
for line in blast:
    if not line.startswith('#'):
        break

#  read in the coordinates
pos = []
for line in blast:
    if line.startswith('#'):
        # stop after the first query
        break

    line = line.rstrip()
    field = line.split()

    # field 6 and 7 are the beginning and end of the query sequence
    field[6] = int(field[6])
    field[7] = int(field[7])
    if field[6] > field[7]:
        # make sure coordinates have smaller value first
        field[6], field[7] = field[7], field[6]

    pos.append({'begin': field[6], 'end': field[7]})

begin = 1
end = 1
for span in sorted(pos, key=lambda p:p['begin']):
    # print('{}\t{}'.format(span['begin'], span['end']))

    if span['begin'] > end:
        # new range - print out missing range, reset begin and end
        print('present\t{}\t{}'.format(begin, end))
        print('missing\t{}\t{}'.format(end+1, span['begin']))
        begin = span['begin']
        end = span['end']

    else:
        # continue current range
        if span['end'] > end:
            end = span['end']

exit(0)
"""=================================================================================================
coverage files from bedtoops genomecov -bga have coverage values at 1 base resolution, making them
very large for plotting.  this script simply filters to values by averaging over a window.  In the
coverage files, regions with the same coverage are given by spans.  Instead of trying to interpolate
an exact window size, i simply stipulate a minimum size over which to accumulate data

example input
tig00000000     0       13      140
tig00000000     13      57      141
tig00000000     57      58      140
tig00000000     58      70      141
tig00000000     70      71      140
tig00000000     71      88      142
tig00000000     88      89      141
tig00000000     89      94      143
tig00000000     94      96      142
tig00000000     96      117     143
tig00000000     117     124     142
id              begin   end     coverage
================================================================================================="""
import sys

minwindow = int(sys.argv[1])
filename = sys.argv[2]

print('file:', filename)
print('window:', minwindow)

try:
    cov = open(filename, 'r')
except:
    print('unable to open coverage file ({})'.format(filename))
    exit(1)

field = []
sum = 0
span = 0
pos_begin = 0;
for line in cov:
    # print(line)
    field  = line.split()
    chr = field[0]
    begin = int(field[1])
    end = int(field[2])
    coverage = int(field[3])

    sum += coverage * (end - begin)
    span += end - begin

    if end - pos_begin > minwindow:
        # write out an average value
        avecov = sum / span
        avepos = (end - pos_begin) / 2
        print('{}\t{}\t{:.2f}\t{:.2f}'.format(pos_begin, end, avepos, avecov))

        pos_begin = end + 1
        span = 0
        sum = 0

# process whatever is left
avecov = sum / span
avepos = (end - pos_begin) / 2
print('{}\t{}\t{:.2f}\t{:.2f}'.format(pos_begin, end, avepos, avecov))



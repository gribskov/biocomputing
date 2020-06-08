"""=================================================================================================
Homework 1

read a fastq file and report the following

total number of lines in the file
number of sequence entries (not the number of lines)
total number of called bases (i.e., total sequence length)

Michael Gribskov     21 January 2019
================================================================================================="""

fastqname = '8044.5k.fastq'

# open the file for reading, terminate with unsuccessful status is file cannot be read
fastq = None        # prevents possibly undefined variable warning
try:
    fastq = open(fastqname, 'r')
except (IOError, OSError):
    print('Error opening file {}'.format(fastqname))
    exit(1)
content = list(fastq)
# read each line and count lines, entries, and bases
n_line = 0
n_entry = 0
n_base = 0
pos = 0
for line in fastq:

    if pos ==1:
        # sequence line
        for base in line:
            try:
                count[base] += 1
            except KeyError:
                count[base] = 1

    pos += 1
    if pos == 4:
        pos = 0

# end of loop over lines of file



# report
print('lines read: {}'.format(n_line))
print('entries read: {}'.format(n_entry))
print('bases read: {}'.format(n_base))

# exit with successful status
exit(0)

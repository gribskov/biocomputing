"""=====================================================================================================================
# HW3_key.py
# Read and store the sequences
# for each sequence, split the sequences after every occurence of 'A' and print the segment in Fasta format, adding an
#   incrementing suffix to indicate the segment of origin
# Repeat with for loop and while loop
# Summarize the total number of sequences found
#
# Michael Gribskov 1/26/2026
====================================================================================================================="""
file = 'scaffolds_short.fa'
try:
    fa = open(file, 'r')
except OSError:
    print(f'Unable to open input file ({file})')

# Read sequences from input file
sequence = {}
for line in fa:
    # read the file line-by-line
    line = line.rstrip()
    if line.startswith('>'):
        # title line, break at first space
        if ' ' in line:
            seqid, junk = line.split(' ', maxsplit=1)
        else:
            seqid = line

        seqid = seqid.replace('>', '')
        sequence[seqid] = ''

    else:
        sequence[seqid] += line

print(f'{len(sequence)} sequences read from {file}')

# print using for loop
total_seg = 0
for s in sequence:
    segment = sequence[s].split('A')
    n_seg = 0
    for seg in segment:
        name = f'{s}_{n_seg}'
        print(f'>{name}\nA{seg}')
        total_seg += 1
        n_seg += 1

print(f'\n{total_seg} subsequences read from {len(sequence)} sequences')

exit(0)

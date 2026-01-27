"""=====================================================================================================================
# HW3_key.py
# Read and store the sequences
# For each sequence, split the sequences after every occurrence of 'A' and print the segment (including the terminal A)
#   in Fasta format, adding an incrementing suffix to indicate the segment of origin
# Repeat with for loop and while loops for the segmentation (inner) loop
# Summarize the total number of sequences and bases found
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
bases_total = 0
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
        bases_total += len(line)

print(f'{len(sequence)} sequences with {bases_total} bases read from {file}\n')

acount = 0
for sid in sequence:
    s = sequence[sid]
    acount += s.count('A')

# Break sequences into segments at each 'A'
total_seg = 0
bases_total = 0
for sid in sequence:
    s = sequence[sid]
    # print(f'sequence: {s}')
    seqlen = len(s)

    begin = 0
    n_seg = 0
    while begin < seqlen:
        end = s.find('A', begin) + 1
        if not end:
            # A not found in remaining segment
            end = seqlen

        name = f'{sid}_{n_seg}'
        # print(f'{name}\t{begin}\t{end}\t{s[begin:end]}')
        print(f'>{name}\n{s[begin:end]}')
        bases_total += end-begin

        n_seg += 1
        total_seg += 1
        begin = end

print(f'\n{total_seg} subsequences with {bases_total} bases read from {len(sequence)} sequences\n')

# Repeat using a for loop for the inner loop
total_seg = 0
bases_total = 0
for sid in sequence:
    s = sequence[sid]
    seqlen = len(s)

    n_pieces = s.count('A')
    if s[-1] != 'A':
        n_pieces += 1
    begin = 0
    n_seg = 0
    for _ in range(n_pieces):
        end = s.find('A', begin) + 1
        if not end:
            # A not found in remaining segment
            end = seqlen

        name = f'{sid}_{n_seg}'
        print(f'>{name}\n{s[begin:end]}')
        bases_total += end-begin

        n_seg += 1
        total_seg += 1
        begin = end

print(f'\n{total_seg} subsequences with {bases_total} bases read from {len(sequence)} sequences\n')

exit(0)

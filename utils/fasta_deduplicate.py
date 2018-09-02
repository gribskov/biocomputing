"""-------------------------------------------------------------------------------------------------
Remove sequences with duplicate names.  Only the first sequence is kept.  Multiple files can be
given for input

Input is two lines, idline and sequence
Output is 100 letters/line

usage
    fasta_deduplicate.py *.fasta > dedup.fasta

2 September 2018    Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import glob
import sys
from sequence.fasta import Fasta

linelen = 100

# default target file name
target = '*.fasta'
if len(sys.argv) > 1:
    target = sys.argv[1]
sys.stderr.write('  target file: {}\n\n'.format(target))

# read all the sequences and store in dictionary
unique_seq = {}
n_file = 0
n_perfile = 0
n_uniqueperfile = 0
n_total = 0
n_unique_total = 0
sys.stderr.write('{}\t{}\t\t{}\n'.format('file', 'per file', 'total'))
for fastafile in glob.glob(target):
    fasta = Fasta()
    fasta.open(fastafile)
    n_file += 1
    n_perfile = 0
    while fasta.next():
        n_perfile += 1
        if fasta.id in unique_seq:
            continue
        else:
            n_uniqueperfile += 1
            unique_seq[fasta.id] = fasta.format(linelen=100)

    n_total += n_perfile
    n_unique_total += n_uniqueperfile
    sys.stderr.write(
        '{}\t{}\t{}\t{}\t{}\t{}\n'.format(n_file, fastafile, n_perfile, n_uniqueperfile, n_total,
                                          n_unique_total))

# write out sequences
for seq in unique_seq:
    sys.stdout.write(unique_seq[seq])

exit(0)

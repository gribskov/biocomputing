"""---------------------------------------------------------------------------------------------------------------------
reformat the output from the apc.pl de-circularization program.
Input is two lines, idline and sequence
Output is 100 letters/line

usage
    fasta_reformat.py *.fasta
---------------------------------------------------------------------------------------------------------------------"""
import glob
import sys
from sequence.fasta import Fasta

linelen = 100

# default target file name
target = '*.fasta'
if len(sys.argv) > 1:
    target = sys.argv[1]
print('  target file:', target)

for fastafile in glob.glob(target):
    # output file
    outfile = fastafile + '.reformatted'
    out = open(outfile, 'w')
    print('  input file:', fastafile, '    output file:', outfile)

    fasta = Fasta()
    fasta.open(fastafile)
    while fasta.next():
        fasta.doc += ' len={}'.format(fasta.length)
        out.write(fasta.format(linelen=linelen))

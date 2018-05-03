"""---------------------------------------------------------------------------------------------------------------------
Remove the Trinity path information from the id line
usage
    fasta_reformat.py *.fasta
---------------------------------------------------------------------------------------------------------------------"""
import glob
import sys
import re
from sequence.fasta import Fasta

linelen = 60

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
        fasta.doc = re.sub(r' path=\[[^]]+\]', '', fasta.doc)
        out.write(fasta.format(linelen=linelen))

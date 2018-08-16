"""-------------------------------------------------------------------------------------------------
convert sequence in raw format to fasta.  Raw format is one sequence/line with no annotation

usage
    raw2fasta <files>

16 August 2018  Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys
import glob
import os
from fasta import Fasta

"""-------------------------------------------------------------------------------------------------
main
-------------------------------------------------------------------------------------------------"""
if __name__ == '__main__':

    inglob = sys.argv[1]
    sys.stderr.write('Input file: {}'.format(inglob))

    for infilename in glob.iglob(inglob):

        try:
            infile = open(infilename, 'r')
        except:
            sys.stderr.write('Unable to open input file ({})\n'.format(infilename))

        base = os.path.basename(infilename)
        base = base.replace('.seq', '')
        sys.stderr.write('Expanded file: {}\t{}\n'.format(infilename, base))
        outfilename = base + '.fasta'
        try:
            outfile = open(outfilename, 'w')
        except:
            yss.stderr.write('Unable to open output file ({})\n'.format(outfilename))

        n = 0
        for seq in infile:
            fasta = Fasta()
            fasta.id = base + '_{}'.format(n)
            fasta.seq = seq.rstrip().upper()
            fasta.doc = 'length={}'.format(fasta.length())
            outfile.write(fasta.format(linelen=100))
            n += 1

        infile.close()
        outfile.close()
        sys.stderr.write('{} sequences in {} written to {}\n'.format(n, infilename, outfilename))

exit(0)

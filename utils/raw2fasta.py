"""-------------------------------------------------------------------------------------------------
convert sequence in raw format to fasta.  Raw format is one sequence/line with no annotation.  Each
file is converted to a fasta formatted file in the local directory.  Sequence characters are
converted to uppercase.  Sequence names are the name of the file with a numeric suffix, e.g. _2

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
    sys.stdout.write('raw2fasta\n')
    sys.stdout.write('\tInput file: {}'.format(inglob))

    for infilename in glob.iglob(inglob):

        # open input fileinfile = None
        try:
            infile = open(infilename, 'r')
        except:
            sys.stderr.write('Unable to open input file ({})\n'.format(infilename))
            exit(1)

        # get base filename and open output file
        base = os.path.basename(infilename)
        base = base.replace('.seq', '')
        sys.stdout.write('\tExpanded file: {}\t{}\n'.format(infilename, base))
        outfilename = base + '.fasta'
        outfile = None
        try:
            outfile = open(outfilename, 'w')
        except:
            sys.stderr.write('Unable to open output file ({})\n'.format(outfilename))
            exit(2)

        # process all sequences in the file
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
        sys.stdout.write(
            '\t{} sequences in {} written to {}\n\n'.format(n, infilename, outfilename))

    # end of loop over files

exit(0)

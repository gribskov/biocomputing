'''
split a sequence in fasta format into chunks
'''
from sequence.fasta import Fasta 
import sys
import argparse
import re

# defaults
maxbases  = 1000000     # maximum number of sequence characters per file
outbase   = "segment"   # default name for output prefix
outsuffix = "fa"        # default suffix (filetype) for output
trim      = ''          # regex string for trimming documentation

commandline = argparse.ArgumentParser( 
    description='Break a fasta file up into segments with no more than a certain'
                'number of sequence characters. Output files are named '
                'PREFIX.n.SUFFIX.'
)
commandline.add_argument( 'fasta_file',
                          help='FastA file to split',
                          type=argparse.FileType('r'))
commandline.add_argument( '--maxbases', 
                          help='maximum number of sequence character per segment', 
                          type=int, 
                          default=str(maxbases))
commandline.add_argument( '--prefix', 
                          help='prefix for output files', 
                          default=outbase )
commandline.add_argument( '--suffix', 
                          help='suffix (filetype) for output files', 
                          default=outsuffix )
commandline.add_argument( '--trim', 
                          help='trim documentation after this regex',
                          default=trim )

# process command line arguments
cl = commandline.parse_args()
maxbases  = cl.maxbases
outbase   = cl.prefix
outsuffix = cl.suffix
outsuffix = outsuffix.lstrip('.')       # remove leading . if present
trim      = cl.trim

print( '\nsplit.py - split fasta file into chunks' )
print( "    fasta file:", cl.fasta_file.name )
print( "    maximum characters:", maxbases )
print( "    output prefix:", outbase )
print( "    output suffix:", outsuffix )
print( "    doc trimmer:", trim )
print( '' )

fasta = Fasta()
fasta.fh = cl.fasta_file

trimre = re.compile( trim )

# initialize counters
base_total   = 0
base_current = 0
n_out        = 0
n_seq        = 0
n_current    = 0

while fasta.next():

	if trimre: fasta.trimDocByRegex( trimre )
	if not n_seq or base_current+fasta.length() > maxbases:
        # if number of bases would be greater than cutoff after adding the new sequence
        # close current output file, open new 
        # report statistics for old file
        # reset current file counters

		try:
			# prevents error on first pass - file is not yet open
			outfile.close()
			print('   ', base_current, 'bases/amino acids', end=' ' )
			print('in', n_current, 'sequences', end=' ' )
			print('written to', outfilename )

		except NameError:
			pass
		
		n_out += 1
		n_current = 0
		base_current = 0

		outfilename = '{0}.{1}.{2}'.format( outbase, n_out, outsuffix )
		outfile = open(outfilename,'w')

	n_seq += 1 
	base_total += fasta.length()

	n_current += 1
	base_current += fasta.length()
	outfile.write( fasta.format() )
	
# report statistics for last file
outfile.close()
print('   ', base_current, 'bases/amino acids', end=' ' )
print('in', n_current, 'sequences', end=' ' )
print('written to', outfilename )

# report overall statistics
print( '\n' )
print( base_total, 'characters from', n_seq, 'sequences written to', n_out, 'files\n' )	

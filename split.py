'''
split a sequence in fasta format into chunks
'''
from sequence.fasta import Fasta 
import sys

# hardwired parameters
outbase = "segment"
outsuffix = "fa"
maxbases = 1000000		# maximum number of sequence characters per file

# input fasta file
datafile = "../../data/card/nucleotide_fasta_protein_homolog_model.fasta"
#datafile = "data/3protein.fa"
print( "fasta file:", datafile )

fasta = Fasta()
fasta.open(datafile)

# initialize counters
base_total   = 0
base_current = 0
n_out        = 0;
n_seq        = 0
n_current    = 0
cutoff       = 0

while fasta.next():
	n_seq += 1 
	base_total += fasta.length()

	if base_total > cutoff:
		try:
			# prevents error on first pass - file is not yet open
			outfile.close()
			print('     ', base_current, 'bases/amino acids in', n_current, 'sequences written to', outfilename)
		except NameError:
			pass
		
		cutoff += maxbases
		n_out += 1
		n_current = 0
		base_current = 0

		outfilename = '{0}.{1}.{2}'.format( outbase, n_out, outsuffix )
		outfile = open(outfilename,'w')

	n_current += 1
	base_current += fasta.length()
	outfile.write( fasta.format() )
	
outfile.close()
print('     ', base_current, 'bases/amino acids in', n_current, 'sequences written to', outfilename)
print( '\n' )
print( base_total, 'bases from', n_seq, 'sequences written to', n_out, 'files' )	
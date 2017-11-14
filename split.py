'''
split a sequence in fasta format into chunks
'''
from sequence.fasta import Fasta 
import sys

#file = "../../data/card/nucleotide_fasta_protein_homolog_model.fasta"
file = "data/1seq.fa"
fasta = Fasta()
#fasta.open(sys.argv)

print( "file:", fasta.filename )
fasta.open(file)

nseq = 0
while fasta.next():
	nseq += 1
	print( "id:", fasta.id )
	print( "doc:", fasta.doc )
	print( "seq:", fasta.seq )
	print( "length:", fasta.length() )
	print( "write:", fasta.format() )

	if nseq > 2: break



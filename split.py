'''
split a sequence in fasta format into chunks
'''
from sequence.fasta import Fasta 
import sys

file = "../../data/card/nucleotide_fasta_protein_homolog_model.fasta"
fasta = Fasta()
#fasta.open(sys.argv)
fasta.open(file)
fasta.read()
print( "id:", fasta.id )
print( "doc:", fasta.doc )
print( "seq:", fasta.seq )

print( "file:", fasta.filename )


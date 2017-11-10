'''
split a sequence in fasta format into chunks
'''
from sequence.fasta import Fasta 
import fileinput as inp
import sys

file = "../../data/card/nucleotide_fasta_protein_homolog_model.fasta"
fasta = Fasta()
#fasta.open(sys.argv)
fasta.open(file)

print( "file:", fasta.filename )


from Bio.Seq import Seq

my_seq = Seq("AGTACACTGGT")
print(my_seq)
print(my_seq.alphabet)

print(my_seq.complement())
print(my_seq.reverse_complement())

print(len(my_seq))
for base in my_seq:
    print(base, end=',')
print()

from Bio.Alphabet import IUPAC

rna = Seq("AGTACACTGGT", IUPAC.unambiguous_rna)
print(rna)
print(rna.reverse_complement())

from Bio import SeqIO
record = SeqIO.read('NC_005816.gb', 'genbank')
print(record)

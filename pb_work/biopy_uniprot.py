# Jilllian Ness [3:46 PM]
import sys
import requests
from Bio import SeqIO
from Bio import ExPASy


params = {"query": "P00750", "format": "fasta", "include": "yes"}
response = requests.get("http://www.uniprot.org/uniprot/", params=params)
# print(response.text)

handle = ExPASy.get_sprot_raw("P00750")
for s in SeqIO.parse(handle, 'swiss'):
    # print(repr(s))
    SeqIO.write(s,sys.stdout,'fasta')




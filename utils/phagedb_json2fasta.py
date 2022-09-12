"""=================================================================================================
phagedb query for all genes
https://phagesdb.org/api/genes/?page=1&page_size=371279

json response looks like
{"count":371279,
 "next":"https://phagesdb.org/api/genes/?page=2",
 "previous":null,
 "results":[
    {"GeneID":"20ES_CDS_23",
     "PhageID":{"PhageID":"20ES",
                "Accession":"KJ410132",
                "Name":"20ES",
                "HostStrain":
                "Mycobacterium",
                "Cluster":"A2"},
     "phams":["56154"],
     "Start":15822,
     "Stop":16230,
     "Length":408,
     "Name":"24",
     "translation":"MTNVFTLDAMREETRKKYQPVKIGLSEDVTVELKPLLKLGKKAREAVADAVKEIEALPDEIDEDDEDSDELMDEVAEKICESIAKVFKLIATSPRKLLAELDTEEEPQIRAELYGAVLRTWMRET QLGEAAPSPN",
     "Orientation":"F",
     "Notes":"b'tail assembly chaperone'"} ...

Michael Gribskov     10 April 2021
================================================================================================="""
import sys
import json
from sequence.fasta import Fasta

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fp = open(sys.argv[1], 'r')
    phage = json.load(fp)

    for gene in phage['results']:
        f = Fasta()
        f.id = gene['GeneID']
        f.seq = gene['translation']
        f.doc = gene['Notes'][2:-1]
        print(f.format(linelen=100))

    exit(0)

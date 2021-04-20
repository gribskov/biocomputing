"""=================================================================================================
Fasta database using MongoDB

Michael Gribskov     20 April 2021
================================================================================================="""
import time
import pymongo
from sequence.fasta import Fasta

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------

mongo = pymongo.MongoClient("mongodb://localhost:27017/")
biocomputing = mongo['biocomputing']
biocomputing.drop_collection('phage')
phage = biocomputing['phage']


fasta = Fasta('C:/Users/michael/Desktop/phage.fa')

start_time = time.perf_counter()
nseq = 0
while fasta.next():
    nseq += 1
    if not nseq %10001:
        break

    result = phage.insert_one( {'_id':fasta.id, 'documentation':fasta.doc, 'sequence':fasta.seq})
    # print(f'\tsequence {fasta.id} inserted as {result.inserted_id}')

end_time = time.perf_counter()
print(f'{phage.count_documents({})} loaded in {end_time-start_time} seconds')

for seq in phage.find():
    print(f'{seq["_id"]}\tlength:{len(seq["sequence"])}')

exit(0)

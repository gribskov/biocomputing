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

fasta_start_time = time.perf_counter()
nseq = 0
all = []
while fasta.next():
    nseq += 1
    # if not nseq % 10001:
    #     break
    all.append({'_id':fasta.id, 'documentation':fasta.doc, 'sequence':fasta.seq})

fasta_end_time = time.perf_counter()
print(f'{nseq} sequences read in {fasta_end_time - fasta_start_time} seconds')

mongo_start_time = time.perf_counter()
result = phage.insert_many(all)
mongo_end_time = time.perf_counter()
print(
    f'{phage.count_documents({})} sequences loaded in {mongo_end_time - mongo_start_time} seconds')
print(f'overall time: {mongo_end_time - fasta_start_time} seconds')

# build index
phage.create_index([('documentation', 'text')])
result = phage.find({'$text':{'$search':'"lysin A"'}}, {'score':{'$meta':'textScore'}})
for seq in result.sort([('score', {'$meta':"textScore"}), ('_id', -1)]):
    print(f'{seq["_id"]} {seq["documentation"]}\t score:{seq["score"]}')

exit(0)

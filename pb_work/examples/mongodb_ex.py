"""=================================================================================================
MongoDB example using pymongo

Michael Gribskov     20 April 2021
================================================================================================="""
import pymongo

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------

mongo = pymongo.MongoClient("mongodb://localhost:27017/")

# database is called biocomputing; collection is called phage
biocomputing = mongo['biocomputing']
phage = biocomputing['phage']

sequence = {'id':'test1', 'documentation':'this is a test', 'sequence':'TCAGAGTCGATGCTGATCTCGC'}
result = phage.insert_one(sequence)

for db in mongo.list_databases():
    print(f'database:{db["name"]}')
    thisdb = mongo[db['name']]
    for col in thisdb.list_collections():
        print(f'\tcollection:{col["name"]}')
        print(f'\tcollection:{col}')

exit(0)

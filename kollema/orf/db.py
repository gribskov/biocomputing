import sys
import json
import pprint
import pymongo
from pymongo import MongoClient
from sequence.trinity import Trinity
from kollema.orf.orf import Orf


class DB:
    """=============================================================================================


    Michael Gribskov     15 May 2021
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.db = MongoClient().peptides

    def load_transcripts(self, fasta, batch=1, clear=True):
        """-----------------------------------------------------------------------------------------
        Load a set of transcripts from a fasta file

        :param fasta: string, path to fasta file
        :return: int, number of transcripts loaded
        -----------------------------------------------------------------------------------------"""
        transcripts = self.db.transcripts
        if clear:
            transcripts.delete_many({})

        trinity = Trinity()
        try:
            trinity.fh = open(fasta, 'r')
        except OSError:
            sys.stderr.write(f'DB:load_transcripts - cannot open file ({fasta})')

        nseq = 0
        datablock = []
        while trinity.next():
            nseq += 1

            trinity.doc = 'len={}'.format(trinity.len)
            trinity.id = trinity.shortid
            # print(trinity.format())
            name = f'{trinity.cluster}_{trinity.component}_{trinity.gene}_{trinity.isoform}'
            orf = Orf()
            orf.sequence = trinity.seq
            nrf = orf.findall(60)
            doc = {'_id':       nseq,
                   'cluster':   f'{trinity.cluster}',
                   'component': f'{trinity.cluster}_{trinity.component}',
                   'gene':      f'{trinity.cluster}_{trinity.component}_{trinity.gene}',
                   'isoform':   name,
                   'sequence':  trinity.seq,
                   'path':      trinity.path,
                   'length':    trinity.len,
                   'rf':        orf.rflist,
                   'nrf':       nrf}

            datablock.append(doc)
            if not nseq % batch:
                transcripts.insert_many(datablock)
                sys.stderr.write(f'loading {nseq}...\n')
                datablock = []

        # load the final block
        transcripts.insert_many(datablock)
        sys.stderr.write(f'loading {nseq}...\n')

        return nseq

    @staticmethod
    def sequences_by_group(collection, level='cluster'):
        """-----------------------------------------------------------------------------------------
        return the set of sequences that match at a specified level:
        cluster, component, or gene (isoform requires no grouping)

        :param level: string, cluster, component, or gene
        :return:
        -----------------------------------------------------------------------------------------"""
        if not level.startswith('$'):
            level = '$' + level

        return collection.aggregate(
            [{"$group":
                  {"_id":     "$cluster",
                   "seqlist": {"$addToSet": "$sequence"},
                   "n":       {"$sum": 1},
                   "name":    {"$addToSet": "$isoform"}
                   }
              }],
            allowDiskUse=True
            )

    def longest_orf_by_gene(self):
        return self.db.transcripts.aggregate([
            {
                # limit for testing
                '$limit': 50
                }, {
                # flatten all reading frames into individual documents, adding the orf sequence
                # and length
                '$unwind': {
                    'path': '$rf'
                    }
                }, {
                '$set': {
                    'rf.seq': {
                        '$substr': [
                            '$sequence', '$rf.begin', {
                                '$subtract': [
                                    '$rf.end', '$rf.begin'
                                    ]
                                }
                            ]
                        },
                    'rf.len': {
                        '$subtract': [
                            '$rf.end', '$rf.begin'
                            ]
                        }
                    }
                }, {
                # group reading frames by gene, get maximum length
                '$group': {
                    '_id':        '$gene',
                    'mxlen':        {
                        '$max': '$rf.len'
                        },
                    'n_isoforms': {
                        '$sum': 1
                        },
                    'iid':        {
                        '$addToSet': '$isoform'
                        },
                    'all':        {
                        '$addToSet': '$rf'
                        }

                    }
                }, {
                # lookup orf with maximum length and report
                '$project': {
                    'length': 1,
                    'mxlen':    1,
                    'mxrf':     {
                        '$arrayElemAt': [
                            '$all', {
                                '$indexOfArray': [
                                    '$all.len', '$mxl'
                                    ]
                                }
                            ]
                        },
                    'mxid':    {
                        '$arrayElemAt': [
                            '$iid', {
                                '$indexOfArray': [
                                    '$all.len', '$mxl'
                                    ]
                                }
                            ]
                        }
                    }
                }, #{
                # '$project': {
                #     '_id': 0
                #     # 'rf':         1,
                #     # 'mxl:': 1,
                #     # 'n_isoforms': 1,
                #     # 'iid':        1,
                #     # 'all':        1
                #     }
                # }
            ])


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    pep = DB()
    print(pep.db.list_collection_names())
    transcripts = pep.db.transcripts
    # # post_id = transcripts.insert_one({'id': 'test'}).inserted_id
    # print(post_id)
    #
    # ntranscript = pep.load_transcripts(
    #     r'A:\mrg\Dropbox\avocado\avocado-R\chile-all\190701_trinity_chile.fasta',
    #     batch=5000,
    #     clear=True)
    # print(f'{ntranscript} transcripts loaded')
    # nn = 0
    # for doc in transcripts.find({}):
    #     nn += 1
    #     print(doc)
    #     if nn > 3:
    #         break

    print('aggregating')
    result = list(pep.longest_orf_by_gene())
    l = list(result)
    print(f'len={len(l)}')
    pp = pprint.PrettyPrinter(indent=2, width=60, compact=False)
    for a in result:
        pp.pprint(a)

    # s = 0
    # for i in range(5):
    #     print(f'i={i}')
    #     s += i

    # result = transcripts.aggregate([
    #     {
    #         '$addFields': {
    #             'length': {
    #                 '$toInt': '$length'
    #                 }
    #             }
    #         }
    #     ])

    # fix length to be int instead of string
    # transcripts.update_many({}, [{'$addFields':{'length':{'$toInt':'$length'}}}])

    # result = transcripts.find({'length': {'$gt': 1500}})
    # i = 0
    # for doc in result:
    #     orf = Orf()
    #     orf.sequence = doc['sequence']
    #     print(doc['isoform'], orf.sequence)
    #     nrf = orf.findall(60)
    #     for rf in sorted(orf.rflist, key=lambda x: x['end'] - x['begin']):
    #         print(rf, rf['end'] - rf['begin'])
    #     transcripts.update_one({'isoform': doc['isoform']},
    #                            {'$set': {'rflist': orf.rflist}})

    # i += 1
    # if i > 10:
    #     break

    # seqs_by_cluster = DB.sequences_by_group(transcripts, level='cluster')
    #
    # for i in seqs_by_cluster:
    #     print(i)

    # db = MongoClient().test
    # db.a.delete_many({})
    # db.a.insert_many([{"a":"dn1"}, {"a":"dn2"}])
    #
    # for i in range(3):
    #     db.a.update_many({"a":"dn1"}, {"$push":{"d":i}})
    # db.a.update_many({"a":"dn2"}, {"$set":{"d":[i for i in range(4)]}})
    # # db.a.update_many({}, {'$push': {'d': [5,6,7]}})
    # # db.a.update_many({}, {'$set': {'e': {'in':7,'out':10}}})
    # # db.a.update_many({}, {'$push': {'f': {'in': 7, 'out': 10}}})
    # result = db.a.find({})
    # pp = pprint.PrettyPrinter(indent=2, width=60, compact=False)
    # print = pp.pprint
    # for doc in result:
    #     # pp.pprint(doc)
    #     print(doc)

    exit(0)

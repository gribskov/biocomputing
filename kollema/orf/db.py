import sys
import json
import pprint
import pymongo
from pymongo import MongoClient
from sequence.trinity import Trinity
from sequence.word import Word
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
            nrf = orf.findall(60, includeseq=True)
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
            # {
            # # limit for testing
            # '$limit': 5000
            # },
            {
                # flatten all reading frames into individual documents, adding the orf sequence
                # and length
                '$unwind': {
                    'path': '$rf'
                    }
                }, {
                '$set': {
                    # 'rf.seq': {
                    #     '$substr': [
                    #         '$sequence', '$rf.begin', {
                    #             '$subtract': [
                    #                 '$rf.end', '$rf.begin'
                    #                 ]
                    #             }
                    #         ]
                    #     },
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
                    'mxlen':      {
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
                    'mxlen':  1,
                    'mxrf':   {
                        '$arrayElemAt': [
                            '$all', {
                                '$indexOfArray': [
                                    '$all.len', '$mxlen'
                                    ]
                                }
                            ]
                        },
                    'mxid':   {
                        '$arrayElemAt': [
                            '$iid', {
                                '$indexOfArray': [
                                    '$all.len', '$mxlen'
                                    ]
                                }
                            ]
                        }
                    }
                },
            {
                '$sort': {'mxlen': -1}
                },
            {
                '$out': 'longorfbygene'
                }
            # {
            # '$project': {
            #     '_id': 0
            #     # 'rf':         1,
            #     # 'mxl:': 1,
            #     # 'n_isoforms': 1,
            #     # 'iid':        1,
            #     # 'all':        1
            #     }
            # }
            ],
            allowDiskUse=True)

    def length_distribution(self, minlen=0, show=True):
        """
        return the orfs greater than minlen and optionally print distribution of length
        :return:
        """
        result = self.db.longorfbygene.find({'mxlen': {'$gte': minlen}})
        if show:
            nval = result.count()
            step = 0.01 * nval
            threshold = 0.0
            n = 0
            for r in result:
                if n >= threshold:
                    frac = f'{n / nval:0.3f}'
                    print(f'{n}\t{frac}\t{r["mxrf"]["len"]}')
                    threshold += step
                n += 1

            frac = f'{n / nval:0.3f}'
            print(f'{n}\t{frac}\t{r["mxrf"]["len"]}')
            print(f'{nval} gene sequences')

        return result


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    load_data = False
    find_long_orfs = False

    pep = DB()
    print(pep.db.list_collection_names())

    if load_data:
        print('loading data')
        transcripts = pep.db.transcripts
        ntranscript = pep.load_transcripts(
            r'A:\mrg\Dropbox\avocado\avocado-R\chile-all\190701_trinity_chile.fasta',
            batch=5000,
            clear=True)
        print(f'{ntranscript} transcripts loaded')

    if find_long_orfs:
        print('finding longest orfs in each gene')
        # transcripts.create_index('gene')
        result = pep.longest_orf_by_gene()
        # result.sort(key=lambda r: r['mxrf']['len'])

    def checkstop(seq):
        nstop = 0
        for pos in range(0,len(seq)-2,3):
            if seq[pos:pos+3] in ('TAA', 'TAG', 'TGA'):
                nstop += 1
        return nstop

    result = pep.length_distribution(1200, show=False)
    coding = Word(size=6)
    coding.init_words(initval=1)
    nword = coding.remove_stops_inframe()
    print(f'{nword} words remaining after removing stops')
    for rf in result:
        check = checkstop(rf["_id"])
        if check > 0:
            print(f'{rf["_id"]}\t{check}')
            print(f'{rf["_id"]}\t{rf["mxlen"]}\t{len(rf["mxrf"]["seq"])}')
        coding.sequence = rf['mxrf']['seq']
        coding.counts_get()

    for word in sorted(coding.count, key=lambda k:coding.count[k], reverse=True):
        # if coding.count[word] == 1:
        #     break

        print(f'{word}\t{coding.count[word]}')
        if word.startswith('TAA') or word.endswith('TAA') or \
            word.startswith('TAG') or word.endswith('TAG') or \
            word.startswith('TGA') or word.endswith('TGA'):
            print("oops")


    # s = 0
    # for i in range(5):
    #     print(f'i={i}')
    #     s += i

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

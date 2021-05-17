import sys
import pymongo
from pymongo import MongoClient
from sequence.trinity import Trinity


class DB:
    """=============================================================================================


    Michael Gribskov     15 May 2021
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.db = MongoClient().peptides

    def load_transcripts(self, fasta, clear=True):
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
        while trinity.next():
            nseq += 1

            trinity.doc = 'len={}'.format(trinity.len)
            trinity.id = trinity.shortid
            # print(trinity.format())
            doc = {'_id':        nseq,
                   'cluster':   trinity.cluster,
                   'component': trinity.component,
                   'gene':      trinity.gene,
                   'isoform':   trinity.isoform,
                   'sequence':  trinity.seq,
                   'path':      trinity.path,
                   'length':    trinity.len}
            post_id = transcripts.insert_one(doc).inserted_id
            # print(post_id)

        return nseq

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    pep = DB()
    print(pep.db.list_collection_names())
    transcripts = pep.db.transcripts
    # post_id = transcripts.insert_one({'id': 'test'}).inserted_id
    # print(post_id)
    #
    # ntranscript = pep.load_transcripts(
    #     r'A:\mrg\Dropbox\avocado\avocado-R\chile-all\190701_trinity_chile.fasta')
    # print(f'{ntranscript} transcripts loaded')


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

    result = transcripts.find({'length':{'$gt':1500}})
    i = 0
    for doc in result:
        print(doc)
        i += 1
        if i>3:
            break

    exit(0)

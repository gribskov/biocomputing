import sys
from math import sqrt, log2
import pprint
import pymongo
from pymongo import MongoClient
from sequence.trinity import Trinity
from sequence.word import Word
from kollema.orf.orf import Orf


class DB:
    """=============================================================================================
    Start MongoDB server:
    C:\Users\michael>..\..\"Program Files"\MongoDB\Server\4.4\bin\mongod --dbpath="a:\MongoDB\Server"

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
        print('finding longest ORFs for each gene')
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

    def frequency_with_orfinfo(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""

        return self.db.longorfbygene.aggregate([
            {
                '$lookup': {
                    'from':         'frequency',
                    'localField':   '_id',
                    'foreignField': '_id',
                    'as':           'freq'
                    }
                }, {
                '$unwind': {
                    'path': '$freq'
                    }
                }, {
                '$project': {
                    'mxlen':      1,
                    'freq.count': 1,
                    'mxrf':       1
                    }
                }
            ])

    def get_long_orfs(self, minlen=0, show=True):

        """-----------------------------------------------------------------------------------------
        return the orfs greater than minlen
        if show is true, print the length distribution

        :param minlen: int, minimum length for ORF
        :param show: boolean, print the distribution
        :return: collection.cursor
        -----------------------------------------------------------------------------------------"""
        # select the longest ORF from each gene equal or longer than minlen
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

    def model_from_counts(self, wordsize=6, initval=1):
        """-----------------------------------------------------------------------------------------
        use the counts in the long ORFs to set up word frequency model
        :param wordsize: int, model word length
        :param initval: int, prior for counts
        :return: Word object with counts
        -----------------------------------------------------------------------------------------"""
        long = self.get_long_orfs(1200, show=False)
        model = Word(size=6)
        model.init_words(initval=1)
        nword = model.remove_stops_inframe()
        print(f'{nword} words remaining after removing stops')
        for rf in long:
            model.sequence = rf['mxrf']['seq']
            model.counts_get()

        return model

    def get_rf_frequencies(self, wordsize=6, background=None, block=10000):
        """-----------------------------------------------------------------------------------------
        use the counts in the long ORFs to set up word frequency model
        :param wordsize: int, model word length
        :param background: Word object, prior for counts, remove unwanted words before using
        :return: Word object with total counts
        -----------------------------------------------------------------------------------------"""
        long = self.get_long_orfs(0, show=False)

        # cloning the background incorporates the prior

        total = background.clone()
        datablock = []
        nseq = 0
        for rf in long:
            rfcount = background.clone()
            nseq += 1
            rfcount.sequence = rf['mxrf']['seq']
            total = rfcount.counts_get()
            rfcount.to_frequency()

            datablock.append({'_id': rf['_id'], 'count': rfcount.count})
            if not nseq % block:
                self.db.frequency.insert_many(datablock)
                datablock = []
                print('.', end='')

        self.db.frequency.insert_many(datablock)

        return total

    def logodds_distribution(self, logodds, wordsize=6, initval=0):
        """-----------------------------------------------------------------------------------------

        :param logodds:
        :param wordsize:
        :param initval:
        :return:
        -----------------------------------------------------------------------------------------"""

        # template for counts from each rf
        template = Word(size=6)
        template.init_words(initval=0)
        nword = template.remove_stops_inframe()

        # calculate logodds score for each rf
        lo_out = open("lo_out.tsv", 'w')
        nrf = 0
        for rf in pep.db.longorfbygene.find({}):
            lo = 0
            rfcount = template.clone()
            rfcount.sequence = rf['mxrf']['seq']
            seqlen = len(rfcount.sequence)
            invsl = rfcount.size / seqlen
            rfcount.count_total = 0
            for word in template.count:
                rfcount.count[word] = invsl
                rfcount.count_total += invsl

            rfcount.counts_get()
            for word in rfcount.count:
                lo += logodds[word] * rfcount.count[word]
            # print(f'{rf["_id"]}\t{lo:.1f}\t{rf["mxrf"]["len"]}')
            lo_out.write(f'{rf["_id"]}\t{lo:.1f}\t{rf["mxrf"]["len"]}\t{lo / rf["mxrf"]["len"]}\n')
            nrf += 1
            if not (nrf % 10000):
                print(nrf)

        lo_out.close()

        return None

    # ==============================================================================================
    # end of DB
    # ==============================================================================================


class Model():
    """=============================================================================================
    centroid models for categories
    ============================================================================================="""

    def __init__(self, template=Word, name=''):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.name = name
        self.current = template.clone()
        self.new = None
        self.len_min = 1000000
        self.len_max = 0
        self.len_ave = 0.0
        self.n = 0
        self.error = 0.0

    def classify_minmax(self, data, frac):
        """---------------------------------------------------------------------------------------------
        score reading frame frequencies vs a single model and classify highest score fraction, frac, as
        model 2, the lowest scoring frac as model 1, and everything else as model 0

        :param data: DB object
        :param frac: float, fraction of top and bottom sequences to classify
        :return: dict, key in rf ID, value is class
        ---------------------------------------------------------------------------------------------"""
        print('assigning to classes')
        model = self.current
        assign = {}
        model_logp = model.to_logfrequency().count
        score = []
        block = 10000
        n = 0

        # calculate the log-likelihood that each ORF belongs to the model
        for freq in data.db.frequency.find({}):
            ll = 0
            count = freq['count']
            for word in freq['count']:
                ll += freq['count'][word] * model_logp[word]
            score.append((freq['_id'], ll))

            n += 1
            if not n % block:
                print(f'+{n}')

        score.sort(key=lambda v: v[1])
        min = int(frac * len(score))
        max = len(score) - min - 1

        assign = {x[0]: 1 for x in score[:min]}  # smallest log likelihood (non-coding)
        assign.update({x[0]: 0 for x in score[min:max]})
        assign.update({x[0]: 2 for x in score[max:]})  # highest log likelihood (coding)

        return assign

    def stat_reset(self):
        """-----------------------------------------------------------------------------------------
        Reset statistics fields in model
        -----------------------------------------------------------------------------------------"""
        self.len_min = 1000000
        self.len_max = 0
        self.len_ave = 0.0
        self.n = 0
        self.error = 0.0

    def stat_addlen(self, rflen):
        """-----------------------------------------------------------------------------------------
        add length of a reading fram to statistics
        -----------------------------------------------------------------------------------------"""
        self.len_min = min(self.len_min, rflen)
        self.len_max = max(self.len_max, rflen)
        self.len_ave += rflen

    def stat_adderror(self, rf_freq):
        """-----------------------------------------------------------------------------------------
        
        :param rf_freq: dict, frequencies of words
        :return:
        -----------------------------------------------------------------------------------------"""
        error = 0
        current = self.current.count
        new = self.new.count
        count_total = self.new.count_total
        count_total = 0
        for word in rf_freq:
            new[word] += rf_freq[word]
            count_total += new[word]
            delta = rf_freq[word] - current[word]
            error += delta * delta

        self.new.count_total = count_total
        error = sqrt(error)
        self.error += error

    def stat_average(self):
        """-----------------------------------------------------------------------------------------
        after all counts have been accumulated, divide error and len_ave by n to get the average.
        divide new frequencies to get an updated model
        :return: 
        -----------------------------------------------------------------------------------------"""
        if self.n > 0:
            self.error /= self.n
            self.len_ave /= self.n
            self.new.count_total = 0
            for word in self.new.count:
                self.new.count[word] /= self.n
                self.new.count_total += self.new.count[word]

    def stat_print(self):
        """-----------------------------------------------------------------------------------------
        format a summary of the model assessment
        :return: 
        -----------------------------------------------------------------------------------------"""
        s = f'{self.name}\t{self.n}\t{self.error:.2g}\t'
        s += f'{self.len_min}\t{self.len_ave:.1f}\t{self.len_max}'

        return s

    # ==============================================================================================
    # end of model
    # ==============================================================================================


def checkstop(seq):
    nstop = 0
    for pos in range(0, len(seq) - 2, 3):
        if seq[pos:pos + 3] in ('TAA', 'TAG', 'TGA'):
            nstop += 1
    return nstop


def composition_from_sequence(dbobj, alphabet='ACGT'):
    """-----------------------------------------------------------------------------------------
    Calculate letter frequencies for all sequences in longest orf by gene
    :param dbobj: Db object
    :param alphabet: string, allowed keys for count
    :return: dict, keys are sequence characters
    -----------------------------------------------------------------------------------------"""
    comp = {i: 0 for i in alphabet}
    comp['total'] = 0
    long = dbobj.get_long_orfs(0, show=False)
    for rf in long:
        for letter in rf['mxrf']['seq']:
            if letter in alphabet:
                comp[letter] += 1
                comp['total'] += 1

    return comp


def update(db, template, model, assign):
    """---------------------------------------------------------------------------------------------
    using category assignments in assign, category 0 is ignored
    Calculate error to current model

    Update model model frequencies
    :param db: DB object
    :param model: list of Model objects
    :param assign: dict of reading frames (keys) and category assignments(values)
    :return:
    ---------------------------------------------------------------------------------------------"""
    new = []
    for m in model:
        m.stat_reset()
        m.new = template.clone()

    orfinfo = db.frequency_with_orfinfo()
    print(f'update stats and categories')
    for orf in orfinfo:
        seq = orf['_id']
        rf_freq = orf['freq']['count']
        m = model[assign[seq]]
        m.stat_addlen(orf['mxlen'])
        m.stat_adderror(rf_freq)
        m.n += 1
        # m.new.add_dict(rf_freq)   # done in adderror

    print(f'\t{model[1].n}\t{model[2].n}')

    for m in model:
        m.stat_average()
        print(m.stat_print())
        m.current = m.new

def classify_nearest(model, data):
    """---------------------------------------------------------------------------------------------
    score reading frame frequencies vs a set of models and reassign to model with maximum
    log-likelihood.  Model zero is reserved for unclassified sequences so it is not tested

    :param frac: list of Model
    :param data: DB object
    :return: dict, key in rf ID, value is class
    ---------------------------------------------------------------------------------------------"""
    print('assigning to nearest class')
    assign = {}
    logp = [None]

    for i in range(1, len(model)):
        # do not include model[0]
        logp.append(model[i].current.to_logfrequency())

    block = 10000
    n = 0

    # calculate the log-likelihood that each ORF belongs to the model
    ccount = [0, 0, 0]
    for freq in data.db.frequency.find({}):
        llmax = -1000000
        llidx = 0
        count = freq['count']
        for i in range(1, len(model)):
            logprob = logp[i].count
            ll = 0
            for word in freq['count']:
                ll += count[word] * logprob[word]

            if ll > llmax:
                llmax = ll
                llidx = i

        assign[freq['_id']] = llidx
        ccount[llidx] += 1

        #todo add counts to new model

        n += 1
        if not n % block:
            print(f'n+{n}\t{ccount[1]}:{ccount[2]}')

    return assign

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # setup template for word counts/frequencies
    # hexamer words minus words with stop codons (which are never present in ORFs)
    count_template = Word(size=6)
    count_template.init_words()
    count_template.remove_stops_inframe()

    # read in transcripts and identify ORFS
    load_data = False
    find_long_orfs = False
    save_frequencies = False

    pep = DB()
    # print(pep.db.list_collection_names())
    # TODO add summary of database collections

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
        rfs = pep.longest_orf_by_gene()
        # result.sort(key=lambda r: r['mxrf']['len'])

    # background distribution from random model, used for prior
    random = count_template.clone()
    composition = composition_from_sequence(pep)
    random.count_random(composition)
    random.to_frequency()

    if save_frequencies:
        # for each long orf, save the hexamer frequencies
        print('calculating word frequencies')
        pep.db.frequency.drop()
        total_freq = pep.get_rf_frequencies(wordsize=6, background=random)

    # set up initial models
    model = [Model(template=count_template.clone(), name='unassigned'),
             Model(template=count_template.clone(), name='noncoding'),
             Model(template=count_template.clone(), name='coding'),
             ]

    # coding model - assume longest ORFs are coding
    model[1].current = pep.model_from_counts(wordsize=6, initval=1)
    model[1].current = model[1].current.to_frequency()

    assign = model[1].classify_minmax(pep, 0.05)
    for i in range(4):
        print(f'\ncycle {i}')
        assign = update(pep, count_template, model, assign)
        classify_nearest(model, pep)

    # logodds = {}
    # for word in coding.count:
    #     logodds[word] = log2(coding_freq[word] / noncoding_freq[word])
    #
    # pep.logodds_distribution(logodds)

    exit(0)

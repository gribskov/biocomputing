"""=================================================================================================
For two blast searches vs the same reference, merge overlapping ranges . This can be used for
to merge multiple transcript or gene predictions

Make sure biocomputing is in your import path

Michael Gribskov     21 March 2025
================================================================================================="""
from blast import Blast
# from ranges.linear_range import Range
from gff.gff import Gff


# import sys
# print(sys.path)

class Feature:
    """=============================================================================================
    A feature is a single linear range it must include
    id: string      id for the range
    begin: int      begin < end, begin is the first position included in the feature (e.g., first
                    base of start codon
    end: int        last base of feature, e.g., last base of stop codon
    source: []      any kind of information to carry along, for instance the original blast hits
                    this allows all the overlapped regions to be identified after merging
    ============================================================================================="""

    def __init__(self, label='', id='', begin=0, end=0):
        self.label = label if label else ''
        self.id = id if id else ''
        self.begin = begin if begin else 0
        self.end = end if end else 0
        self.source = []


class Range:
    """=============================================================================================
    Generalization from the old range class. A range is a list of dictionaries. The dictionaries
    must include the following required fields:
            id: string
            begin: int
            end: int
    All other fields are essentially passive payload
    ============================================================================================="""

    def __init__(self, sortfunc=None):
        """-----------------------------------------------------------------------------------------
        sort        sorting function to use for detecting overlaps

        -----------------------------------------------------------------------------------------"""
        self.features = []

        if not sortfunc:
            self.sort = self.sort_default

    @staticmethod
    def sort_default(data):
        """-----------------------------------------------------------------------------------------
        Default sort is first by id, then by begin

        :param data: dict       dictionary to sort
        -----------------------------------------------------------------------------------------"""
        return (data.id, data.begin)

    def overlap(self, mindist=0):
        """-----------------------------------------------------------------------------------------
        Merge overlapping ranges if the distance between the end of the current range and the
        beginning of the next range is < mindist. After merging, the features list may contain
        elements with the value None (the ranges that were merged)

        :param mindist: int     ranges with a gap > mindist do not overlap (positive integer)
        :return: int            number of ranges
        -----------------------------------------------------------------------------------------"""
        if len(self.features) < 2:
            return

        self.features.sort(key=self.sort)

        feature_n = 1
        current = self.features[0]
        remove = []

        for next in self.features[1:]:
            if current.id != next.id:
                # different reference sequences, overlap is impossible, next becomes current
                current = next
                continue

            if next.begin - current.end - 1 <= mindist:
                # merge ranges, current stays the same
                current.begin = min(current.begin, next.begin)
                current.end = max(current.end, next.end)
                for s in next.source:
                    current.source.append(s)
                # merged ranges are added to remove list, can't remove here because it changes the list being
                # iterated
                remove.append(next)
            else:
                # no overlap, next becomes current
                current = next
                feature_n += 1

        if remove:
            # remove merged eatures
            for merged in remove:
                self.features.remove(merged)

        return feature_n

    def merge(self, other):
        """-----------------------------------------------------------------------------------------
        Add the features in other to self

        :param other: Range     an existing Range object
        :return: int            number of features in self after merge
        -----------------------------------------------------------------------------------------"""
        for f in other.features:
            self.features.append(f)

        return len(self.features)


def read_and_filter_blast(infile, columns, evalue=1e-5, pid=95):
    """---------------------------------------------------------------------------------------------
    to make the forward and reverse regions not overlap, matches in the opposite direction have 'r'
    appended to the name of the reference sequence name. This is appropriate for stranded RNASeq
    libraries where there can be different transcripts on each strand.

    :param infile:
    :param evalue:
    :param pid:
    :return:
    ---------------------------------------------------------------------------------------------"""
    # lrange will be a list of all features for this query
    all_ranges = Range()

    i = 0
    for block in blast_query_set(infile, columns):
        # block is the blast result lines for a single query
        lrange = Range()
        for value in block:
            if value['pid'] < pid or value['evalue'] > evalue:
                # both pid and evalue must be true or we go on the next entry in the block
                continue

            # passed both pid and evalue thresholds
            thisrange = Feature(id=value['sid'], begin=value['sbegin'], end=value['send'])
            thisrange.source.append(value)

            if value['reverse']:
                # mark sequences on reverse strand and reverse begin and end so begin < end
                # reverse=True is reverse strand
                thisrange.id += 'r'
                thisrange.begin, thisrange.end = thisrange.end, thisrange.begin

            lrange.features.append(thisrange)

        # all features for this query are in lrange, check for overlaps and merge
        # save merged features in all_ranges
        lrange.overlap(mindist=10000)
        all_ranges.merge(lrange)
        i += 1

    return all_ranges


def blast_query_set(infile, columns):
    """---------------------------------------------------------------------------------------------
    Generator that reads and returns a block of results with the same query sequence. This function
    sets reverse=True if send < sbegin

    :param infile:
    :param columns:
    :return:
    ---------------------------------------------------------------------------------------------"""
    msplit = len(columns) - 1
    cname = list(columns.keys())

    blast = open(infile, 'r')

    block = None
    qid_old = None
    for line in blast:
        if line.startswith('#'):
            continue

        fields = line.rstrip().split(maxsplit=msplit)
        values = {cname[i]: fields[i] for i in range(len(columns))}
        for c in values:
            # change strings to numbers where needed
            if columns[c] == 'i':
                values[c] = int(values[c])
            elif columns[c] == 'f':
                values[c] = float(values[c])

        values['reverse'] = False  # forward strand
        if values['send'] < values['sbegin']:
            values['reverse'] = True  # reverse strand

        if values['qid'] != qid_old:
            if block:
                # yield sorted(block, key=lambda b: (b['reverse'], b['sid'], min(b['sbegin'], b['send'])))
                yield block
            qid_old = values['qid']
            block = [values]
            continue

        block.append(values)

    if block:
        # yield sorted(block, key=lambda b: (b['reverse'], b['sid'], min(b['sbegin'], b['send'])))
        yield block

    blast.close()
    return


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    trinityfile = 'data/c16c31.trinity.stuberosum.blastn'

    searchfields = {'qid':   's', 'qlen': 'i', 'qbegin': 'i', 'qend': 'i',
                    'sid':   's', 'sbegin': 'i', 'send': 'i',
                    'qcov':  'i', 'allen': 'i', 'pid': 'f',
                    'score': 'f', 'evalue': 'f', 'stitle': 's'}

    # read in blastfiles and and filter by evalue
    # hits may be exons if the search is transcripts vs genome so aggregate hits within maxsep
    # on the same subject sequence into a single range
    trinity = read_and_filter_blast(trinityfile, searchfields, evalue=1e-5, pid=95)

    exit(0)

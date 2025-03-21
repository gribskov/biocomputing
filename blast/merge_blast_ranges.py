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

class Range:
    """=============================================================================================
    Try to generalize from the old range class. A range is a list of dictionaries. The dictionaries
    must include the following required fields:
            id: string
            begin: int
            end: int
    All other fields are essentially passive payload
    ============================================================================================="""

    def __init__(self, sortfunc= None):
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
        return (data['id'], data['begin'])


def read_and_filter_blast(infile, columns, evalue=1e-5, pid=95):
    """---------------------------------------------------------------------------------------------

    :param infile:
    :param evalue:
    :param pid:
    :return:
    ---------------------------------------------------------------------------------------------"""
    lrange = Range()

    for block in blast_query_set(infile, columns):
        if len(block) > 5:
            print('stop')
        for values in block:
            if values['pid'] < pid or values['evalue'] > evalue:
                continue

        thisrange = {'id': values['sid'], 'begin': values['sbegin'], 'end': values['send'], 'source': [values]}
        if values['reverse']:
            # reverse=True is reverse strand
            thisrange['id'] += 'r'
            thisrange['begin'], thisrange['end'] = thisrange['end'], thisrange['begin']

        lrange.features.append(thisrange)

    # TODO add sorting here
    return


def blast_query_set(infile, columns):
    """---------------------------------------------------------------------------------------------
    Generator that reads and returns a block of results with the same query sequence

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

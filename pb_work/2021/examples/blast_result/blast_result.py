"""=================================================================================================
Manipulate a tabular output Blast result. Hits are store internally as a list of dictionaries.

Michael Gribskov     15 February 2021
================================================================================================="""
import sys


class BlastResult:

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------
        BlastResult constructor
        :param filename: str, readable blast tabular output (format 6 or 7)
        -----------------------------------------------------------------------------------------"""
        self.filename = filename
        self.fh = None
        self.hits = []

        if filename:
            self.opensafe()

    def opensafe(self):
        """-----------------------------------------------------------------------------------------
        Open a file for reading with error trapping.  the filename should be in the BlastResult
        object. Exit status 1 if file open fails

        Usage
            BlastResult.opensafe()

        :return: filehandle, filehandle of opened file
        -----------------------------------------------------------------------------------------"""
        try:
            self.fh = open(self.filename, 'r')
        except (IOError, OSError):
            sys.stderr.write('BlastResult::opensafe - error opening file {}'.format(self.filename))
            exit(1)

        return self.fh

    def read_and_parse(self):
        """-----------------------------------------------------------------------------------------
        Read and parse the entire blast result.  Store in BlastResult.hits as a list of dictionaries
        Dictionary keys are:
            query_id
            subject_id
            percent_id
            query_coverage
            bit_score
            evalue

        :return: int, number of hits
        -----------------------------------------------------------------------------------------"""
        nhits = 0
        for line in self.fh:
            nhits += 1
            field = line.rstrip().split()
            self.hits.append({'qseqid':    field[0],
                              'sseqid':    field[1],
                              'pident':    float(field[2]),
                              'length':    int(field[3]),
                              'mismatch':  int(field[4]),
                              'gapopen':   int(field[5]),
                              'qstart':    int(field[6]),
                              'qend':      int(field[7]),
                              'sstart':    int(field[8]),
                              'ssend':     int(field[9]),
                              'bit_score': float(field[10]),
                              'evalue':    float(field[11]),
                              })
        return nhits


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    blastx = BlastResult('genome.small_uniref50.dmndblastx')
    nhits = blastx.read_and_parse()
    sys.stderr.write('{} blast hits read from {}\n'.format(nhits, blastx.filename))

    # for hit in sorted(blastx.hits, key=lambda h: h['sseq:id']):
    #     sys.stdout.write('{}\n'.format(hit))
    #
    # subject = {}
    # subj_idx = {}
    # n_subject = 0
    # n_hit = 0
    # for hit in blastx.hits:
    #     if hit['subject_id'] in subject:
    #         # a known subject
    #         subject[hit['subject_id']] += 1
    #     else:
    #         # previously unknown subject
    #         subject[hit['subject_id']] = 1
    #         n_subject += 1
    #         subj_idx[hit['subject_id']] = []
    #
    #     subj_idx[hit['subject_id']].append(hit)
    #     n_hit += 1
    #
    # sys.stderr.write('{} unique subjects found\n'.format(n_subject))
    #
    # for s in sorted(subj_idx, key=lambda s: (-len(subj_idx[s]), s)):
    #     # sys.stderr.write('{}: {}\n'.format(s, subj_idx[s]))
    #     sys.stderr.write('{}\n'.format(s))
    #     for h in subj_idx[s]:
    #         sys.stderr.write('\t{}\n'.format(h))
    #
    # print('{} => {}'.format(id(blastx.hits[0]), blastx.hits[0]))
    # s = 'C1MXJ7_MICPC'
    # hit = subj_idx[s]
    # print('{} => {}'.format(id(hit[0]), hit[0]))

    exit(0)

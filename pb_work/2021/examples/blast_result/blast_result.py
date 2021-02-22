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
            self.hits.append({'qseqid':   field[0],
                              'sseqid':   field[1],
                              'pident':   float(field[2]),
                              'length':   int(field[3]),
                              'mismatch': int(field[4]),
                              'gapopen':  int(field[5]),
                              'qstart':   int(field[6]),
                              'qend':     int(field[7]),
                              'sstart':   int(field[8]),
                              'ssend':    int(field[9]),
                              'evalue':   float(field[10]),
                              'bitscore': float(field[11]),
                              })
        return nhits


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    blastx = BlastResult('genome.small_uniref50.dmndblastx')
    nhits = blastx.read_and_parse()
    sys.stderr.write('{} blast hits read from {}\n'.format(nhits, blastx.filename))

    # for hit in sorted(blastx.hits, key=lambda h: h['sseqid']):
    #     sys.stdout.write('{}\t\t{}\t\t{}\n'.format(hit['qseqid'], hit['sseqid'], hit['evalue']))

    # count unique sseqid and number of hits for each, store the index of each hit for each sseqid
    subject = {}
    subj_idx = {}
    n_subject = 0
    n_hit = 0
    for hit in blastx.hits:
        if hit['sseqid'] in subject:
            # a known subject
            subject[hit['sseqid']] += 1
        else:
            # previously unknown subject
            subject[hit['sseqid']] = 1
            n_subject += 1
            subj_idx[hit['sseqid']] = []

        subj_idx[hit['sseqid']].append(n_hit)
        n_hit += 1

    sys.stderr.write('{} unique subjects found\n'.format(n_subject))
    nn = 0
    for s in sorted(subj_idx, key=lambda s: (-len(subj_idx[s]), s)):
        sys.stderr.write('{}:{}\n'.format(s, len(subj_idx[s])))
        for h in subj_idx[s]:
            sys.stderr.write('\t{}\t{}\t{}\t{}\n'.
                             format(h, blastx.hits[h]['sseqid'], blastx.hits[h]['qseqid'],
                                    blastx.hits[h]['evalue']))
        nn += 1
        if nn > 5:
            break

    exit(0)

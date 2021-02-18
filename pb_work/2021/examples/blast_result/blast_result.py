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
        Read and parse the entire blast result.  Store in BlastResult.hits as a lsit of dictionaries
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
            self.hits.append({'query_id':       field[0],
                              'subject_id':     field[1],
                              'percent_id':     float(field[2]),
                              'query_coverage': float(field[3]),
                              'bit_score':      int(field[4]),
                              'evalue':         float(field[5]),
                              }
                             )
        return nhits


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    blastx = BlastResult('../../../../blast/data/diamond.blastx')
    nhits = blastx.read_and_parse()
    sys.stderr.write('{} blast hits read from {}\n'.format(nhits, blastx.filename))

    for hit in sorted(blastx.hits, key=lambda h: h['subject_id']):
        sys.stdout.write('{}\n'.format(hit))

    exit(0)

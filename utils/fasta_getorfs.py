"""=================================================================================================
fasta_getorfs.py
Get protein open reading frames from nucleic acid files in fasta format
Usually used to get protein sequences from predicted transcripts

Michael Gribskov    13 December 2018
================================================================================================="""
import sys
# import re
import argparse
from sequence.fasta import Fasta


class Orf():
    """=============================================================================================
    Open reading frame object.

    ============================================================================================="""

    def __init__(self, transcript=None):
        """-----------------------------------------------------------------------------------------
        Orf constructor

        :param transcript, Fasta object, DNA sequence for reading frame extraction
        -----------------------------------------------------------------------------------------"""
        self.transcript = transcript
        self.orf = []  # list of ORFs
        self.min_len = 50  # minimum length for orfs

    def get(self, direction='fr', frame=(0, 1, 2)):
        """-----------------------------------------------------------------------------------------
        Add ORFs >= min_len to ORF list (self.orf)
        TODO: correctly calculate position for forward and reverse
        TODO: change direction to +- instead of rf

        :param direction, string - use 'f' for forward, 'r' for reverse
        :param frame, tuple with a list of the reading frames to translate
        :return: n_added, integer
        -----------------------------------------------------------------------------------------"""
        n_orf = 0

        for dir in direction:

            for f in frame:
                trans = self.transcript.translate(frame=f, direction=dir)
                begin = f

                for pep in trans.seq.split('*'):
                    peplen = len(pep)
                    end = begin + 3 * peplen

                    if peplen > args.minlen:
                        n_orf += 1

                        # print('\t',pep)

                        this_orf = {}
                        this_orf['id'] = '{}_{}'.format(self.transcript.id, n_orf)
                        this_orf['direction'] = dir
                        this_orf['frame'] = f
                        this_orf['begin'] = begin
                        this_orf['end'] = end
                        this_orf['length'] = peplen
                        this_orf['sequence'] = pep

                        self.orf.append(this_orf)

                    begin = end + 3

        return len(self.orf)

    def write_as_fasta(self, fh, n=None):
        """-----------------------------------------------------------------------------------------
        Write to a file in fasta format, if n is defined, write only the specified ORF in the list

        :param fh, open filehandle for writing
        :param n:
        :return:
        -----------------------------------------------------------------------------------------"""
        fasta = Fasta()
        for orf in self.orf:
            fasta.id = orf['id']
            fasta.doc = 'strand={} frame={} begin={} end={}'.format(orf['direction'], orf['frame'], orf['begin'], orf['end'])
            fasta.seq = orf['sequence']
            fh.write(fasta.format(linelen=60))
            fh.write('\n')

        return

    def write_as_tabular(self, n=None):
        """-----------------------------------------------------------------------------------------
        Write to a file in tabular format, if n is defined, write only the specified ORF in the list

        :param n:
        :return:
        -----------------------------------------------------------------------------------------"""
        pass

        return


def arguments_get():
    minlen_default = 50
    cl = argparse.ArgumentParser(description='Get ORFs from transcript sequences')
    cl.add_argument('-m', '--minlen', type=int, default=minlen_default,
                    help='minimum length open reading frame to report')
    cl.add_argument('fasta_in', type=argparse.FileType('r'))

    return cl.parse_args()


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = arguments_get()

    fasta = Fasta(fh=args.fasta_in)
    fasta.read()
    print(fasta.format())

    orf = Orf(fasta)
    orf.get()

    out = open('a.a', 'w')
    orf.write_as_fasta(out)

    sys.stderr.write('done\n')

exit(0)

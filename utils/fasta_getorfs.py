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


class Orf:
    """=============================================================================================
    Open reading frame object.

    As usual with sequence features, numbering can be confusing.  The sequence strings are treated
    as normal zero-based python strings.  When the begin/end coordinates are written out, they are
    changed to the conventional one-based convention you biologists are more familiar with.

    Ranges begin at the first base of the first codon, and end at the last base of a non-stop codon.
    The stop codon is not included.

    Synopsis

    fasta = Fasta(fh=args.fasta_in)
    nsequence = 0
    npeptide = 0
    npeptide_total = 0
    while fasta.next():
        nsequence += 1

        orf = Orf(fasta)
        orf.min_len = args.minlen
        orf.get()

        npeptide = orf.write_as_fasta(sys.stdout)
        npeptide_total += npeptide
        orf.write_as_tabular(sys.stderr)
        sys.stderr.write('\n{} peptide sequences written from {}\n'.format(npeptide, fasta.id))

    ============================================================================================="""

    def __init__(self, transcript=None):
        """-----------------------------------------------------------------------------------------
        Orf constructor

        :param transcript, Fasta object, DNA sequence for reading frame extraction
        -----------------------------------------------------------------------------------------"""
        self.transcript = transcript
        self.orf = []  # list of ORFs
        self.min_len = 50  # minimum length for orfs

    def get(self, direction='+-', frame=(0, 1, 2)):
        """-----------------------------------------------------------------------------------------
        Add ORFs >= min_len to ORF list (self.orf)
        All ORFs are numbered, even if they are shorter than the minimum length.  This ensures that 
        the naming will be the same, even if run with different parameters

        coordinates are calculated in terms of the forward strand with the first base as position 1
        this is the same convention as used for GFF/GTF

        :param direction, string - use '+' for forward, '-' for reverse
        :param frame, tuple with a list of the reading frames to translate
        :return: n_added, integer
        -----------------------------------------------------------------------------------------"""
        n_orf = 0
        seqlen = len(self.transcript.seq)

        for strand in direction:

            for f in frame:
                trans = self.transcript.translate(frame=f, direction=strand)
                begin = f
                # len_eff = len(trans.seq) * 3
                # if strand == '-':
                #     begin = (len(trans.seq) - f) % 3

                for pep in trans.seq.split('*'):
                    n_orf += 1
                    peplen = len(pep)
                    end = begin + 3 * peplen

                    if peplen >= self.min_len:
                        start = begin + 1
                        stop = end
                        if strand == '-':
                            start = seqlen - end + 1
                            stop = seqlen - begin

                        this_orf = {'id': '{}_{}'.format(self.transcript.id, n_orf),
                                    'direction': strand,
                                    'frame': f,
                                    'begin': start,
                                    'end': stop,
                                    'length': peplen,
                                    'sequence': pep}

                        self.orf.append(this_orf)

                    begin = end + 3
                    # end of loop of ORFS
                    # end of loop over frames
                    # end of loop over strands (directions)

        return len(self.orf)

    def write_as_fasta(self, fh, n=None):
        """-----------------------------------------------------------------------------------------
        Write to a file in fasta format, if n is defined, write only the specified ORF in the list

        :param fh, open filehandle for writing
        :param n: integer, index of ORF to write (not implemented)
        :return: n
        -----------------------------------------------------------------------------------------"""
        fasta = Fasta()
        nwritten = 0

        if n:
            orf = self.orf[n]
            fasta.id = orf['id']
            fasta.doc = 'len={} strand={} frame={} begin={} end={}'. \
                format(orf['length'], orf['direction'], orf['frame'], orf['begin'], orf['end'])
            fasta.seq = orf['sequence']
            fh.write(fasta.format(linelen=60))
            fh.write('\n')
            nwritten = 1

        for orf in self.orf:
            fasta.id = orf['id']
            fasta.doc = 'len={} strand={} frame={} begin={} end={}'. \
                format(orf['length'], orf['direction'], orf['frame'], orf['begin'], orf['end'])
            fasta.seq = orf['sequence']
            fh.write(fasta.format(linelen=60))
            fh.write('\n')
            nwritten += 1

        return nwritten

    def write_as_tabular(self, fh, n=None):
        """-----------------------------------------------------------------------------------------
        Write to a file in tabular format, if n is defined, write only the specified ORF in the list

        :param fh, open filehandle for writing
        :param n:  integer, index of ORF to write (not implemented)
        :return: nwritten number writtedn
        -----------------------------------------------------------------------------------------"""
        nwritten = 0
        if n:
            orf = self.orf[n]
            fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                orf['id'], orf['direction'], orf['frame'], orf['begin'],
                orf['end'], orf['length'], orf['sequence']))
            nwritten = 1

        else:
            for orf in self.orf:
                fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    orf['id'], orf['direction'], orf['frame'], orf['begin'],
                    orf['end'], orf['length'], orf['sequence']))
                nwritten += 1

        return nwritten

    # ====== end of class Orf ======================================================================


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
    peptide_maxlen = 20000
    args = arguments_get()
    sys.stderr.write('\nfasta_getorfs - Get ORFs from transcript sequences\n')
    sys.stderr.write('\tinput transcript file: {}\n'.format(args.fasta_in.name))
    sys.stderr.write('\tminimum ORF length: {}\n\n'.format(args.minlen))

    fasta = Fasta(fh=args.fasta_in)
    nsequence = 0
    npeptide = 0
    npeptide_total = 0

    lenhist = [0 for _ in range(peptide_maxlen)]
    lentotal = 0

    longhist = [0 for _ in range(peptide_maxlen)]
    longtotal = 0

    while fasta.next():
        nsequence += 1

        orf = Orf(fasta)
        orf.min_len = args.minlen
        orf.get()

        longest = True
        for pep in sorted(orf.orf,key=lambda k: k['length'],reverse=True):
            if longest:
                longhist[pep['length']] += 1
                longtotal += 1
                longest = False

            lenhist[pep['length']] += 1
            lentotal += 1


        # npeptide = orf.write_as_fasta(sys.stdout)
        npeptide_total += npeptide
        # orf.write_as_tabular(sys.stderr)
        sys.stderr.write('\n{} peptide sequences written from {}\n'.format(npeptide, fasta.id))

    sys.stderr.write('\n{} transcripts read from {}\n'.format(nsequence, args.fasta_in.name))
    sys.stderr.write('{} peptides with len > {} extracted\n'.format(npeptide_total, args.minlen))

    cumulative = 0
    print('\nall ORFs\n{:>8s}{:>8s}{:>10s}{:>8s}{:>10s}'.format('# len', 'n', 'frac', 'sum', 'sumfrac'))
    for i in range(peptide_maxlen):
        if lenhist[i]:
            cumulative += lenhist[i]
            print('{:8d}{:8d}{:10.4g}{:8d}{:10.4f}'.format(i,
                                                           lenhist[i], lenhist[i] / lentotal,
                                                           cumulative, cumulative / lentotal))

    cumulative = 0
    print('\nlongest ORFs\n{:>8s}{:>8s}{:>10s}{:>8s}{:>10s}'.format('# len', 'n', 'frac', 'sum', 'sumfrac'))
    for i in range(peptide_maxlen):
        if longhist[i]:
            cumulative += longhist[i]
            print('{:8d}{:8d}{:10.4g}{:8d}{:10.4f}'.format(i,
                                                           longhist[i], longhist[i] / longtotal,
                                                           cumulative, cumulative / longtotal))

exit(0)

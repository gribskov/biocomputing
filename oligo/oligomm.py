"""=================================================================================================
get a count of kmers in a sequence.  The real goal is to identify the infrequent kmers and use them
to probabilistically generate longer kmers unlikely to occur in the sequence.

TODO: write oligos to file
TODO: check for duplicates?
TODO: bias towards more extreme AT/GC content?
================================================================================================="""
import sys
from math import log10 as log
import random
import itertools


class Kmer:
    """=============================================================================================
    Generates and holds the kmer counts and probabilities

    ============================================================================================="""

    def __init__(self, k=8):
        """-----------------------------------------------------------------------------------------
        kmer constructor
        :param k: length of words
        -----------------------------------------------------------------------------------------"""
        self.k = k
        self.alphabet = 'ACGT'
        self.count = {}
        self.p = {}
        self.w = {}
        self.setupWords()
        self.total = 0
        self.pmin = 1.0
        self.pmax = 0.0

    def reset(self):
        """-----------------------------------------------------------------------------------------
        Reset the kmer counts, probabilities, and number of kmers
        :return: True
        -----------------------------------------------------------------------------------------"""
        self.list = []
        self.count = {}
        self.p = {}
        self.w = {}
        self.pmin = 1.0
        self.pmax = 0.0

        return True

    def setupWords(self, prior=1):
        """-----------------------------------------------------------------------------------------
        initialize all of the keys and give a count equal to prior
        :param prior:
        :return: nwords
        -----------------------------------------------------------------------------------------"""
        n = 0
        self.list = []
        for word in itertools.product(self.alphabet, repeat=self.k):
            kstr = ''.join(word)
            self.list.append(kstr)
            self.count[kstr] = prior
            self.p[kstr] = 0.0
            self.w[kstr] = 1.0
            n += 1

        return n

    def fromFasta(self, fasta):
        """-----------------------------------------------------------------------------------------
        count kmers in a Fasta file.  Some words may be excluded because they contain non-alphabet
        characters.

        :param file: open file handle for fasta file
        :return: number of kmers counted
        -----------------------------------------------------------------------------------------"""

        # count words in sequence
        self.total = 0
        nseq = 0
        seq = ''  # for overlap
        trace = ''
        for line in fasta:
            line = line.rstrip()
            if line.startswith('>'):
                nseq += 1
                trace = nseq
                # print('\nn', nseq)
                seq = ''  # no overlap between sequences
            else:
                line = seq + line
                for i in range(len(line) - self.k + 1):
                    try:
                        self.count[line[i:i + self.k]] += 1

                        # screen trace: TODO make interval an option
                        if not self.total % 10000000:
                            print()
                            print('{:13d} '.format(self.total), end='')
                        if not self.total % 100000:
                            print('{}'.format(trace), end='')
                            trace = '.'
                        self.total += 1

                    except KeyError:
                        # kmers with non ACGT letters
                        # sys.stderr.write('unknown kmer: {}\n'.format(line[i:i + KMER]))
                        pass
                # restore for debugging
                # if self.total > 10000000:
                #     break

                sys.stdout.flush()
                seq = line[-self.k:-1]  # this is the KMER-1 letters that overlap the next line

        return self.total
        # end of fromFasta

    def updateProb(self):
        """-----------------------------------------------------------------------------------------
        Update the probabilities from the current counts.  Assumes the count in self.total is
        correct.

        :return: pmin, pmax
        -----------------------------------------------------------------------------------------"""
        self.pmin = 1.0
        self.pmax = 0.0

        for word in self.count:
            self.p[word] = self.count[word] / self.total
            self.pmin = min(self.pmin, self.p[word])
            self.pmax = max(self.pmax, self.p[word])

        return self.pmin, self.pmax

    def tableWrite(self, kmerout):
        """-----------------------------------------------------------------------------------------
        write the kmer table to a file
        :param kmerout: file opened by argparse
        :return: kmers written
        -----------------------------------------------------------------------------------------"""
        total = 0
        for word in self.count:
            kmerout.write('{}\t{}\t{:.4g}\n'.format(word, self.count[word], self.p[word]))
            total += self.count[word]

        return total

    def tableRead(self, kmerin, reset=True):
        """-----------------------------------------------------------------------------------------
        Read in a kmer table written by tableWrite.  If reset is True, the kmer count and
        probability are reset
        :param file: file opened by argparse
        :return: integer, number of words read
        -----------------------------------------------------------------------------------------"""

        # reset the object data if requested
        if reset:
            self.reset()

        # read the file
        n = 0
        for line in kmerin:
            if line.startswith(('#', '!')):
                continue

            kmer, count, p = line.split()
            self.count[kmer] = int(count)
            pf = float(p)
            self.p[kmer] = pf
            self.pmin = min(pf, self.pmin)
            self.pmax = max(pf, self.pmax)
            if n == 0:
                self.k = len(kmer)
            n += 1

        self.total = n
        return n


# end of Kmer class ================================================================================
import sys
import argparse


def commandLine(default):
    """---------------------------------------------------------------------------------------------
    Get command line options

    :return:
    ---------------------------------------------------------------------------------------------"""
    commandline = argparse.ArgumentParser(
        description='calculate infrequent oligos based on kmer frequencies in an input'
                    'Fasta file. Oligos are written Fasta format to stdout.'

    )
    commandline.add_argument('--fasta',
                             help='FastA file to split',
                             type=argparse.FileType('r'),
                             )

    commandline.add_argument('--kmer',
                             help='maximum number of sequence character per segment',
                             type=int,
                             default=str(default['kmer'])
                             )

    commandline.add_argument('--noligo',
                             help='number of oligos to generate',
                             type=int,
                             default=str(default['noligo'])
                             )

    commandline.add_argument('--table',
                             help='precalculated kmer table',
                             type=argparse.FileType('r'),
                             )

    commandline.add_argument('--koutput',
                             help='output kmer table',
                             type=argparse.FileType('w')
                             )

    return commandline.parse_args()


def printParameters(clargs):
    """---------------------------------------------------------------------------------------------
    Print a report on command line parameters to stderr.

    :param clargs: command line arguments from argparse.parse_args()
    :return: True
    ---------------------------------------------------------------------------------------------"""
    sys.stderr.write('\n')
    if clargs.fasta or clargs.table:
        for k in clargs.__dict__:
            if k.startswith('__') or clargs.__dict__[k] == None:
                continue
            if k == 'fasta' or k == 'table':
                sys.stderr.write('{}: {}\n'.format(k, clargs.__dict__[k].name))
            else:
                sys.stderr.write('{}: {}\n'.format(k, clargs.__dict__[k]))

    else:
        sys.stderr.write('Input error: fasta (--fasta) table or kmer table (--table) must be provided\n')
        exit(1)

        sys.stderr.write('\n')
    return True


def fractionGC(seq, alphabet):
    """---------------------------------------------------------------------------------------------
    Calculate fraction GC in a sequence. non ACGT letters are ignored

    :param seq:
    :return: fraction GC
    ---------------------------------------------------------------------------------------------"""
    count = {x: 0 for x in alphabet}
    all = 0
    for base in seq:
        count[base] += 1

        all += 1

    return (count['G'] + count['C']) / all


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    default = {'kmer': 4, 'noligo': 1000}
    cl = commandLine(default)
    printParameters(cl)

    kmer = Kmer(cl.kmer)

    nwords = 0
    infile = ''
    if cl.table:
        # read kmers from file
        nwords = kmer.tableRead(cl.table)
        print('{} {}'.format(kmer.k, kmer.total))
        infile = cl.table.name
    else:
        # calculate kmers from fasta
        # convert word counts to probabilities and write out
        nwords = kmer.fromFasta(cl.fasta)
        pmin, pmax = kmer.updateProb()
        infile = cl.fasta.name

    sys.stderr.write('\ntotal {}mer words read from {}: {}\n'.format(kmer.k, infile, kmer.total))
    sys.stderr.write('\n    pmin pmax: {:10.3g}\t{:10.3g}\n'.format(kmer.pmin, kmer.pmax))
    sys.stderr.write('log pmin pmax: {:10.3g}\t{:10.3g}\n'.format(log(kmer.pmin), log(kmer.pmax)))

    if cl.output:
        nwritten = kmer.tableWrite(cl.output)
        sys.stderr.write('\n{} kmers written to {}'.format(nwritten, cl.output))

    # weight counts - ad hoc weighting function
    # w = bias**(logP - logPmax)
    # the larger the bias, the more the weights favor infrequent words

    wmin = 1.0
    wmax = 0.0
    wsum = 0.0
    bias = 0.5
    cutoff = log(kmer.pmax)
    threshold = {}
    # for sorted list
    #  for word in sorted(kmer.count, key=lambda k: kmer.count[k]):
    for word in kmer.list:
        # TODO it seems like the following should be kmer.p not kmer.count
        # TODO need to rethink the whole calculation but need longer kmer data
        kmer.w[word] = bias ** (cutoff - log(kmer.p[word]))
        wsum += kmer.w[word]
        wmin = min(wmin, kmer.w[word])
        wmax = max(wmax, kmer.w[word])
        threshold[word] = wsum
        # print('{:6d} {:10}{:10.3g}'.format(nw, word, kmer.count[word]))

    # print weighted range
    sys.stderr.write('w: {:10.3g}{:10.3g}\n'.format(wmin, wmax))

    # could overlap words by KMER-1, but this is unnecessary.  try overlap = 3
    overlap = 3
    omer = {}
    for word in itertools.product(kmer.alphabet, repeat=overlap):
        # generate all overlap len words
        oword = ''.join(word)
        omer[oword] = {'list': [], 't': [], 'w': 0.0}

    # calculate the transitions between overlapping words: 'list' is the list of overlapping full
    # words, 't' has the threshold value for each word (cumulative weigth), and 'w' has the total
    # weight summed over the overlapping words
    for word in kmer.list:
        lap = word[:overlap]
        omer[lap]['list'].append(word)
        omer[lap]['w'] += kmer.w[word]
        omer[lap]['t'].append(omer[lap]['w'])

    # construct an oligo

    for n in range(cl.noligo):
        # select a weighted start point
        r = random.random() * wsum
        found = ''
        for word in kmer.list:
            found = word
            if r < threshold[word]:
                break

        oligo = found
        weight = kmer.w[found]

        nsteps = 8
        for step in range(nsteps):
            # select an overlapping next word
            lap = oligo[-overlap:]
            r = random.random() * omer[lap]['w']
            found = ''
            for i in range(len(omer[lap]['list'])):
                found = omer[lap]['list'][i]
                if r < omer[lap]['t'][i]:
                    break
            oligo += found[-overlap:]

        gc = fractionGC(oligo, kmer.alphabet)
        sys.stdout.write('>n{} GC={:.3}\n'.format(n, gc))
        sys.stdout.write(oligo)

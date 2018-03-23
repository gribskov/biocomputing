"""=================================================================================================
get a count of kmers in a sequence.  The real goal is to identify the infrequent kmers and use them
to probabilistically generate longer kmers unlikely to occur in the sequence.

TODO: read saved kmer distribution
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
        self.kmer = {}
        self.total = 0
        self.p = {}
        self.pmin = 1.0
        self.pmax = 0.0

    def setupWords(self, prior=1):
        """-----------------------------------------------------------------------------------------
        initialize all of the keys and give a count equal to prior
        :param prior:
        :return: nwords
        -----------------------------------------------------------------------------------------"""
        n = 0
        for word in itertools.product(self.alphabet, repeat=self.k):
            kstr = ''.join(word)
            self.kmer[kstr] = prior
            self.p[kstr] = 0.0
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
                        self.kmer[line[i:i + self.k]] += 1

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

        for word in self.kmer:
            self.p[word] = self.kmer[word] / self.total
            self.pmin = min(self.pmin, self.p[word])
            self.pmax = max(self.pmax, self.p[word])

        return self.pmin, self.pmax

    def tableWrite(self, file):
        """-----------------------------------------------------------------------------------------
        write the kmer table to a file
        :param file: file name
        :return: kmers written
        -----------------------------------------------------------------------------------------"""
        kmerout = None
        try:
            kmerout = open(file, 'w')
        except OSError:
            sys.stderr.write('kmer.tablewrite - unable to open output file ({})'.format(file))

        total = 0
        for word in self.kmer:
            kmerout.write('{}\t{}\t{:.4g}\n'.format(word, self.kmer[word], self.p[word]))
            total += self.kmer[word]

        return total


# end of kKmer class ===============================================================================
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

    commandline.add_argument('--output',
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
        sys.stderr('Input error: fasta (--fasta) table or kmer table (--table) must be provided\n')
        exit(1)

        sys.stderr.write('\n')
    return True


def fractionGC(seq):
    """---------------------------------------------------------------------------------------------
    Calculate fraction GC in a sequence. non ACGT letters are ignored

    :param seq:
    :return: fraction GC
    ---------------------------------------------------------------------------------------------"""
    count = {}
    all = 0
    for base in seq:
        try:
            count[base] += 1
        except KeyError:
            count[base] = 1

        all += 1

    return (count['G'] + count['C']) / all


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    default = {'kmer': 8, 'noligo': 1000}
    cl = commandLine(default)
    printParameters(cl)

    kmer = Kmer(cl.kmer)
    nwords = kmer.setupWords()
    print('{} words initialized'.format(nwords))

    nwords = kmer.fromFasta(cl.fasta)
    print('\ntotal {}mer words read from {}: {}'.format(cl.kmer, sys.argv[1], nwords))

    # convert word counts to probabilities and write out
    pmin, pmax = kmer.updateProb()
    print('\n    pmin pmax: {:10.3g}\t{:10.3g}'.format(pmin, pmax))
    print('log pmin pmax: {:10.3g}\t{:10.3g}'.format(log(pmin), log(pmax)))

    nwritten = kmer.tableWrite(sys.argv[2])
    print('\n{} kmers written to {}'.format(nwritten, sys.argv[2]))

    # weight counts - ad hoc weighting function
    # w = bias**(logP - logPmax)
    # the larger the bias, the more the weights favor infrequent words

    wmin = 1.0
    wmax = 0.0
    wsum = 0.0
    cutoff = log(pmax)
    threshold = {}
    # for sorted list
    #  for word in sorted(kmer.kmer, key=lambda k: kmer.kmer[k]):
    for word in kmer.kmer:
        kmer.kmer[word] = 4 ** (cutoff - log(kmer.kmer[word]))
        wsum += kmer.kmer[word]
        wmin = min(wmin, kmer.kmer[word])
        wmax = max(wmax, kmer.kmer[word])
        threshold[word] = wsum
        # print('{:6d} {:10}{:10.3g}'.format(nw, word, kmer.kmer[word]))

    # print weighted range
    print('w: {:10.3g}{:10.3g}'.format(wmin, wmax))

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
    for word in kmer.kmer:
        lap = word[:overlap]
        omer[lap]['list'].append(word)
        omer[lap]['w'] += kmer.kmer[word]
        omer[lap]['t'].append(omer[lap]['w'])

    # construct an oligo

    for n in range(0, 1000):
        # select a weighted start point
        r = random.random() * wsum
        found = ''
        for word in kmer.kmer:
            found = word
            if r < threshold[word]:
                break

        oligo = found
        weight = kmer.kmer[found]

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

        gc = fractionGC(oligo)
        print('>n{} GC={:.3}'.format(n, gc))
        print(oligo)

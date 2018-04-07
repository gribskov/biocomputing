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
        self.wcum = {}
        self.setupWords()
        self.total = 0
        self.pmin = 1.0
        self.pmax = 0.0
        self.wmin = 1.0
        self.wmax = 0.0

    def reset(self):
        """-----------------------------------------------------------------------------------------
        Reset the kmer counts, probabilities, and number of kmers

        :return: True
        -----------------------------------------------------------------------------------------"""
        self.list = []
        self.count = {}
        self.p = {}
        self.w = {}
        self.wcum = {}

        self.total = 0
        self.pmin = 1.0
        self.pmax = 0.0
        self.wmin = 1.0
        self.wmax = 0.0
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
                trace = '\n{}'.format(nseq)
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

        # check if probabilities and weights are available and write either 1, 2 or 3 values
        total = 0
        if self.pmax > 0 and self.wmax > 0:
            # write count, p, weight
            for word in self.count:
                kmerout.write('{}\t{}\t{:.4g}\t{:.4g}\n'.format(
                    word, self.count[word], self.p[word], self.w[word]))
                total += self.count[word]

        elif self.pmax > 0:
            # write count, p
            for word in self.count:
                kmerout.write('{}\t{}\t{:.4g}\n'.format(word, self.count[word], self.p[word]))
                total += self.count[word]

        else:
            # $ write count only
            for word in self.count:
                kmerout.write('{}\t{}\n'.format(word, self.count[word]))
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

            field = line.split()
            k = field[0]
            self.count[k] = int(field[1])
            self.p[k] = float(field[2])
            self.pmin = min(field[2], self.pmin)
            self.pmax = max(field[2], self.pmax)

            if len(field) > 3:
                self.w[k] = float(field[3])
                self.wmin = min(field[3], self.wmin)
                self.wmax = max(field[3], self.wmax)

            n += 1

        # assume k is the length of the first word
        keys = self.count.keys()
        self.k = len(keys[0])

        self.total = n
        return n

    def weightNegExp(self, exp=1.0):
        """-----------------------------------------------------------------------------------------
        Similar to dirichlet wighting, inverted.  values range [0 - 1]
        weight = 1 - f(kmer) ** exp

        :parameter exp: exponent to which the base is taken
        :return: float, min, max weight
        -----------------------------------------------------------------------------------------"""
        wmax = 0
        wmin = 1
        for k in self.p:
            # w = 1.0 - self.p[k] ** exp
            w = self.p[k] ** exp
            wmax = max(w, wmax)
            wmin = min(w, wmin)
            self.w[k] = w
            print('{:.4g}\t{}'.format(w, k))

        self.wmin = wmin
        self.wmax = wmax

        return wmin, wmax

    def weightExp(self, exp=1.0):
        """-----------------------------------------------------------------------------------------
        Similar to dirichlet wighting, inverted.  values range [0 - 1]
        weight = f(kmer) ** exp

        :parameter exp: exponent to which the base is taken
        :return: float, min, max weight
        -----------------------------------------------------------------------------------------"""
        wmax = 0
        wmin = 1
        for k in self.p:
            # w = 1.0 - self.p[k] ** exp
            w = self.p[k] ** exp
            wmax = max(w, wmax)
            wmin = min(w, wmin)
            self.w[k] = w

        self.wmin = wmin
        self.wmax = wmax

        return wmin, wmax

    def weightNormalize(self):
        """-----------------------------------------------------------------------------------------
        Divide weights by sum.  update wmin and wmax and cumulative weights, self.wcum
        :return: float, min and max weight
        -----------------------------------------------------------------------------------------"""
        wsum = 0.0
        for k in self.w:
            wsum += self.w[k]

        wmin = 1.0
        wmax = 0.0
        wcum = 0.0
        for k in self.list:
            w = self.w[k] / wsum
            wmax = max(w, wmax)
            wmin = min(w, wmin)
            self.w[k] = w

            wcum += w
            self.wcum[k] = wcum

        self.wmin = wmin
        self.wmax = wmax

        return wmax, wmin

    def randomByWeight(self):
        """-----------------------------------------------------------------------------------------
        select a ranom kmer using the current stored weights.  Assumes weights have been normalized
        to [0.0, 1.0].  See weightNormalize

        :return: string random kmer, float weight
        -----------------------------------------------------------------------------------------"""
        r = random.random()
        found = ''
        for word in self.list:
            found = word
            if r > self.wcum[word]:
                break

        return found, self.w[found]




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

    commandline.add_argument('--overlap',
                             help='overlap between kmer words',
                             type=int,
                             default=str(default['overlap'])
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

    default = {'kmer': 8, 'noligo': 1000, 'overlap':3}
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

    print('\ntotal {}mer words read from {}: {}'.format(kmer.k, infile, kmer.total))

    # weight counts - ad hoc weighting function
    # w = bias**(logP - logPmax)
    # the larger the bias, the more the weights favor infrequent words, save this for future ref

    kmer.weightExp(-2.0)
    kmer.weightNormalize()
    print('\n    pmin pmax: {:10.3g}\t{:10.3g}'.format(kmer.pmin, kmer.pmax))
    print('log pmin pmax: {:10.3g}\t{:10.3g}'.format(log(kmer.pmin), log(kmer.pmax)))
    print('\n    wmin wmax: {:10.3g}\t{:10.3g}'.format(kmer.wmin, kmer.wmax))

    if cl.output:
        nwritten = kmer.tableWrite(cl.output)
        print('\n{} kmers written to {}'.format(nwritten, cl.output))

    omer = Kmer(cl.overlap)
    # calculate the transitions between overlapping words: the overlap transition is simply a list
    # of shorter kmers
    # TODO not right.  for each overlap word create a kmer object with just the overlapping kmers.
    for word in kmer.count:
        lap = word[:overlap]
        omer.count[lap] += kmer.count[word]
        omer.p[lap] += kmer.w[word]
        omer.w[lap] += kmer.w[word]

    # construct oligos

    for n in range(cl.noligo):
        # select a weighted start point
        oligo, weight = kmer.randomByWeight()

        nsteps = 8
        for step in range(nsteps):
            # select an overlapping next word
            lap = oligo[-overlap:]
            r = random.random() * omer.w[lap]
            found = ''
            for i in range(len(omer[lap]['list'])):
                found = omer[lap]['list'][i]
                if r < omer[lap]['t'][i]:
                    break
            oligo += found[-overlap:]

        gc = fractionGC(oligo, kmer.alphabet)
        print('>n{} GC={:.3}'.format(n, gc))
        print(oligo)

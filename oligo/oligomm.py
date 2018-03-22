"""=================================================================================================
get a count of kmers in a sequence.  The real goal is to identify the infrequent kmers and use them
to probabilistically generate longer kmers unlikely to occur in the sequence.
================================================================================================="""
import sys
from math import log10 as log
import random
import itertools


class Kmer:
    """=============================================================================================
    Holds the kmer distribution

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

    def fromFasta(self, file):
        """-----------------------------------------------------------------------------------------
        count kmers in a Fasta file.  Some words may be excluded because they contain non-alphabet
        characters.

        :param file:
        :return: number of kmers counted
        -----------------------------------------------------------------------------------------"""
        fasta = None
        try:
            fasta = open(file, 'r')
        except OSError:
            sys.stderr.write('kmer.fromFasta - unable to open output file ({})'.format(file))

        # count words in sequence
        self.total = 0
        nseq = 0
        seq = ''  # for overlap
        for line in fasta:
            line = line.rstrip()
            if line.startswith('>'):
                nseq += 1
                print('\nn', nseq)
                seq = ''  # no overlap between sequences
            else:
                line = seq + line
                for i in range(len(line) - KMER + 1):
                    try:
                        self.kmer[line[i:i + KMER]] += 1
                        if not self.total % 10000000:
                            print()
                            print('{:13d} '.format(self.total), end='')
                        if not self.total % 100000:
                            print('.', end='')
                        self.total += 1

                    except KeyError:
                        # kmers with non ACGT letters
                        # sys.stderr.write('unknown kmer: {}\n'.format(line[i:i + KMER]))
                        pass
                # TODO: remove for production
                if self.total > 10000000:
                    break

                sys.stdout.flush()
                seq = line[-KMER:-1]    # this is the KMER-1 letters that overlap the next line

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


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    KMER = 8
    kmer = Kmer(KMER)
    nwords = kmer.setupWords()
    print('{} words initialized'.format(nwords))

    nwords = kmer.fromFasta(sys.argv[1])
    print('\ntotal {}mer words read from {}: {}'.format(KMER, sys.argv[1], nwords))

    # convert word counts to probabilities and write out
    pmin, pmax = kmer.updateProb()
    print('    pmin pmax: {:10.3g}\t{:10.3g}'.format(pmin, pmax))
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
    print('w: {:10.3g}{:10.3g}'.format(pmin, pmax))

    overlap = 3
    omer = {}
    for word in itertools.product(kmer.alphabet, repeat=overlap):
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

"""=================================================================================================
get a count of kmers in a sequence.  The real goal is to identify the infrequent kmers and use them
to probabilistically generate longer kmers unlikely to occur in the sequence.
================================================================================================="""
import sys
import pprint
from math import log2 as log
import random
import itertools
from sequence.fasta import Fasta


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

    def setupWords(self, prior=1):
        """-----------------------------------------------------------------------------------------
        initialize all of the keys and give a count equal to prior
        :param prior:
        :return: nwords
        -----------------------------------------------------------------------------------------"""
        n = 0
        for word in itertools.product(self.alphabet, repeat=self.k):
            self.kmer[''.join(word)] = prior
            n += 1

        return n


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    sys.stdout.flush()

    fasta = Fasta()
    fasta.open(sys.argv[1])

    KMER = 6
    kmer = Kmer(KMER)
    nwords = kmer.setupWords()
    print('{} words initialized'.format(nwords))

    # count words in sequence
    total = 0
    nseq = 0
    while fasta.next():
        nseq += 1
        if not nseq % 1000:
            print('.', end='', flush=True)
        if not nseq % 20000:
            print()

        for i in range(len(fasta.seq) - KMER + 1):
            kmer.kmer[fasta.seq[i:i + KMER]] += 1
            total += 1
        if nseq > 100:
            break

    print('\ntotal {}mer words: {}'.format(KMER, total))
    pp = pprint.PrettyPrinter(indent=4)
    # pp.pprint(kmer.kmer)

    # convert word counts to probabilities
    pmin = 1.0
    pmax = 0.0
    for word in kmer.kmer:
        kmer.kmer[word] /= total
        pmin = min(pmin, kmer.kmer[word])
        pmax = max(pmax, kmer.kmer[word])

    print('p: {:10.3g}{:10.3g}'.format(pmin, pmax))
    print('p: {:10.3g}{:10.3g}'.format(log(pmin), log(pmax)))

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
        kmer.kmer[word] = 4**(cutoff - log(kmer.kmer[word]))
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
        omer[oword] = {'list': [], 't':[], 'w': 0.0}

    # calculate the transitions between overlapping words: 'list' is the list of overlapping full
    # words, 't' has the threshold value for each word (cumulative weigth), and 'w' has the total
    # weight summed over the overlapping words
    for word in kmer.kmer:
        lap = word[:overlap]
        omer[lap]['list'].append(word)
        omer[lap]['w'] += kmer.kmer[word]
        omer[lap]['t'].append(omer[lap]['w'])

    # construct an oligo

    #select a weighted start point
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



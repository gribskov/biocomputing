""""=================================================================================================
count of kmers in a sequence. Calculates probabilities and probability weighted kmer sampling

TODO: add synopsis/testing
TODO: check for duplicates?
TODO: bias towards more extreme AT/GC content?

2 March 2018     Michael Gribskov
================================================================================================="""
import sys
import random
import itertools


class Kmer():
    """=============================================================================================
    Generates and holds the kmer counts and probabilities

    ============================================================================================="""

    def __init__(self, k=8, setup=True):
        """-----------------------------------------------------------------------------------------
        kmer constructor.
        :param k: length of words
        -----------------------------------------------------------------------------------------"""
        self.k = k
        self.alphabet = 'ACGT'
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

        if setup:
            self.setupWords()

    def reset(self):
        """-----------------------------------------------------------------------------------------
        Reset the kmer counts, probabilities, and number of kmers. k and alphabet are not reset.

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

    def add(self, attr, key, value=1):
        """-----------------------------------------------------------------------------------------
        safe add function.  if try fails, dict element is created
        :param attr: object attribute name
        :return: True
        -----------------------------------------------------------------------------------------"""
        try:
            self.__dict__[attr][key] += value
        except KeyError:
            self.__dict__[attr][key] = value

        return True

    def fromFasta(self, fasta, interval=1000000):
        """-----------------------------------------------------------------------------------------
        count kmers in a Fasta file.  Some words may be excluded because they contain non-alphabet
        characters.

        :param fasta: open file handle for fasta file
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
                # sys.stderr.write('\nn', nseq)
                seq = ''  # no overlap between sequences
            else:
                line = seq + line
                for i in range(len(line) - self.k + 1):
                    try:
                        self.count[line[i:i + self.k]] += 1

                        # screen trace:
                        if not self.total % interval:
                            sys.stderr.write('\n{:13d} '.format(self.total))
                            if not self.total % interval / 10:
                                sys.stderr.write('{}'.format(trace))
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
            for word in self.list:
                kmerout.write('{}\t{}\t{:.4g}\t{:.4g}\n'.format(
                    word, self.count[word], self.p[word], self.w[word]))
                total += self.count[word]

        elif self.pmax > 0:
            # write count, p
            for word in self.list:
                kmerout.write('{}\t{}\t{:.4g}\n'.format(word, self.count[word], self.p[word]))
                total += self.count[word]

        else:
            # $ write count only
            for word in self.list:
                kmerout.write('{}\t{}\n'.format(word, self.count[word]))
                total += self.count[word]

        return total

    def tableRead(self, kmerin, reset=True):
        """-----------------------------------------------------------------------------------------
        Read in a kmer table written by tableWrite.  If reset is True, the kmer count and
        probability are reset
        :param kmerin: file opened by argparse
        :param reset: if True, clear attributes before loading
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
            self.list.append(k)
            self.count[k] = int(field[1])
            self.p[k] = float(field[2])
            self.pmin = min(self.p[k], self.pmin)
            self.pmax = max(self.p[k], self.pmax)

            if len(field) > 3:
                self.w[k] = float(field[3])
                self.wmin = min(self.w[k], self.wmin)
                self.wmax = max(self.w[k], self.wmax)

            n += 1

        # assume k is the length of the first word
        self.k = len(list(self.count)[0])

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
        for k in self.list:
            # w = 1.0 - self.p[k] ** exp
            w = self.p[k] ** exp
            wmax = max(w, wmax)
            wmin = min(w, wmin)
            self.w[k] = w
            # sys.stderr.write('{:.4g}\t{}'.format(w, k))

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
        for k in self.list:
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
        Divide weights by sum.  update self.wmin and self.wmax and cumulative weights, self.wcum
        :return: float, min and max weight
        -----------------------------------------------------------------------------------------"""
        wsum = 0.0
        for k in self.list:
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
        select a random kmer using the current stored weights.  Assumes weights have been normalized
        to [0.0, 1.0].  See weightNormalize

        :return: string random kmer, float weight
        -----------------------------------------------------------------------------------------"""
        r = random.random()
        found = ''
        for word in self.list:
            found = word
            if r < self.wcum[word]:
                break

        return found, self.w[found]


# end of Kmer class ================================================================================

# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    exit(0)

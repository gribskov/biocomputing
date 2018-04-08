"""=================================================================================================
Example of EM for promoter identification

Michael Gribskov     08 April 2018
================================================================================================="""
import sys
import random


class pssm():
    """=============================================================================================
    position specific scoring matrix class.
    ============================================================================================="""

    def __init__(self, size=6, seqfile=sys.stdin):
        """-----------------------------------------------------------------------------------------
        Pssm constructor
        :param seqfile: string, name of the file containing the sequences
        -----------------------------------------------------------------------------------------"""
        self.size = size
        self.pssm = []
        self.alphabet = ''
        self.composition = {}
        self.sequence = []
        self.expected = []
        self.likelihood = []

        self.sequenceRead(seqfile)
        self.compositionCalc()

    def initRandom(self, p):
        """-----------------------------------------------------------------------------------------
        Fill pssm with background composition +/- random value [-p, +p)
        :param p: fraction of random noise to add to composition
        :return: float min, float max
        -----------------------------------------------------------------------------------------"""
        # conversion for random [0.0, 1.0) to [-p, +p)
        # p = = (r-0.5)2p = 2pr - p
        twop = p * 2
        # initialize pssm and column sums
        colsum = [0 for k in range(self.size)]
        self.pssm = [[0 for k in range(self.size)] for a in range(len(self.alphabet))]

        #  calculate column sums as we go
        i = 0
        colsum = [0 for k in range(self.size)]
        for letter in self.alphabet:
            # rows are sequence characters
            for pos in range(self.size):
                r = random.random() * twop - p
                self.pssm[i][pos] = self.composition[letter] + r
                colsum[pos] += self.composition[letter] + r
            i += 1

        # convert pssm values to probabilities by dividing by column sums
        valmin = 1.0
        valmax = 0.0
        i = 0
        for letter in self.alphabet:
            # rows are sequence characters
            for pos in range(self.size):
                self.pssm[i][pos] /= colsum[pos]
                valmin = min(self.pssm[i][pos], valmin)
                valmax = max(self.pssm[i][pos], valmax)
            i += 1

        return valmin, valmax

    def expectation(self):
        """-----------------------------------------------------------------------------------------
        Expectation step
        Calculate that each sequence position is a site. scores for each sequence are not normalized
        to sum to 1 across the sequences
        :return: True
        -----------------------------------------------------------------------------------------"""
        seqnum = 0
        self.expected = []
        for s in self.sequence:
            exp = []
            self.expected.append(exp)
            for pos in range(len(s) - self.size + 1):
                i = self.alphabet.find(s[pos])
                p = self.pssm[i][0]
                for w in range(1,self.size):
                    i = self.alphabet.find(s[pos+w])
                    p *= self.pssm[i][w]
                exp.append(p)
            seqnum += 1

        return True

    def ml(self):
        """-----------------------------------------------------------------------------------------
        Maximization step
        Calclulate the maximum likelihood pattern based on the probabilities that each sequence
        position is a site.  The ML estimate is simply the probability weighted occurrence of the
        letters in each position of the pssm
        :return: True
        -----------------------------------------------------------------------------------------"""
        # reinitialize pssm
        for a in range(len(self.alphabet)):
            for w in range(self.size):
                self.pssm[a][w] = 0.0

        colsum = [0 for k in range(self.size)]
        seqnum = 0
        for s in self.sequence:
            for pos in range(len(s) - self.size + 1):
                for w in range(0,self.size):
                    i = self.alphabet.find(s[pos+w])
                    self.pssm[i][w] += self.expected[seqnum][pos+w]
                    colsum[w] += self.pssm[i][w]
            seqnum += 1

        # renormalize columns by dividing by colsums.  pssm should be probability  afterwards

        for a in range(len(self.alphabet)):
            for w in range(self.size):
                self.pssm[a][w] /= colsum[w]

        return True

    def table(self):
        """-----------------------------------------------------------------------------------------
        Return a list of lists with the current pssm

        :return: list, 2D list of sequence letters (rows) and positons(columns)
        -----------------------------------------------------------------------------------------"""
        pass

    def sequenceRead(self, seqfile):
        """-----------------------------------------------------------------------------------------
        read a set of sequences, one per line with no non-sequence characters
        duplicate seqeunces are not included
        no checking for non-sequence characters

        -----------------------------------------------------------------------------------------"""
        try:
            seq = open(seqfile, 'r')
        except:
            sys.stderr.write(
                'pssm.sequenceRead - unable to open sequence file ({})'.format(seqfile))

        self.nseq = 0
        for line in seq:
            line = line.rstrip()
            if line not in self.sequence:
                self.sequence.append(line)
                self.nseq += 1

        seq.close()
        return self.nseq

    def compositionCalc(self):
        """-----------------------------------------------------------------------------------------
        Calculate fractional composition of sequences.  sets composition and alphabet
        :return: number of letters
        -----------------------------------------------------------------------------------------"""
        count = 0
        for s in self.sequence:
            for letter in s:
                try:
                    self.composition[letter] += 1
                except KeyError:
                    self.composition[letter] = 1
                count += 1

        for letter in self.composition:
            self.composition[letter] /= count

        self.alphabet = ''.join(list(self.composition))
        return count

    def convergence(self, table):
        """-----------------------------------------------------------------------------------------
        Calculate the mean sum of squared differences between the current pssm and table

        :return: float, mean squared difference
        -----------------------------------------------------------------------------------------"""
        mse = 1.0

        return mse

    # End of pssm class ============================================================================


# ==================================================================================================
#
# ==================================================================================================
if __name__ == '__main__':
    # read promoter sequences and composition
    model = pssm(6, seqfile=sys.argv[1])

    # create random pssm from random frequencies +/- random
    p_rand = 0.05
    model.initRandom(p_rand)

    # EM - until convergence
    model.expectation()
    model.ml()
    pssm_current = model.table
    converged = False
    target = 0.01
    while not converged:
        # calculate position probabilities for each sequence (expectation)
        model.expectation()

        # calculate maximimum likelihood pssm
        model.ml()

        mse = model.convergence(pssm_current)
        if mse < target:
            converged = True

    # report probability distributions

    # report pssm

    exit(0)

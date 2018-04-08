"""=================================================================================================
Example of EM for promoter identification

Michael Gribskov     08 April 2018
================================================================================================="""
import sys

class pssm():
    """=============================================================================================
    position specific scoring matrix class.
    ============================================================================================="""

    def __init__(self, size=6, seqfile=sys.stdin):
        """-----------------------------------------------------------------------------------------
        Pssm constructor
        :param seqfile: string, name of the file containing the sequences
        -----------------------------------------------------------------------------------------"""
        pass

    def initRandom(self, p):
        """-----------------------------------------------------------------------------------------

        :param p: fraction of random noise to add to composition
        :return: float min, float max
        -----------------------------------------------------------------------------------------"""

    def expectation(self):
        """-----------------------------------------------------------------------------------------
        Expectation step
        Calculate that each sequenc position is a site

        -----------------------------------------------------------------------------------------"""
        pass

    def ml(self):
        """-----------------------------------------------------------------------------------------
        Maximization step
        Calclulate the maximum likelihood pattern based on the probabilities that each sequence
        position is a site.  The ML estimate is simply the probability weighted occurrence of the
        letters in each position of the pssm

        -----------------------------------------------------------------------------------------"""
        pass

    def table(self):
        """-----------------------------------------------------------------------------------------
        Return a list of lists with the current pssm

        :return: list, 2D list of sequence letters (rows) and positons(columns)
        -----------------------------------------------------------------------------------------"""
        pass

    def sequenceRead(self, seqfile):
        """-----------------------------------------------------------------------------------------
        read a set of sequences, one per line with no non-sequence characters

        -----------------------------------------------------------------------------------------"""
        pass

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
    model = pssm(6, seqfile=sys.argv[0])

    # create random pssm from random frequencies +/- random
    p_rand = 0.1
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

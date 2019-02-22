import sys


# ==================================================================================================
class Fasta:

    def __init__(self, filename):
        pass

    def next(self):
        return True

    def composition(self):
        return True

    def mw(self):
        return True


# end of class Fasta ===============================================================================

def mw_average(mwlist):
    """---------------------------------------------------------------------------------------------
    Calculate the average molecular weight. Returns 0.0 if the list is empty.

    :param mwlist: list of float, mw of each protein
    :return: float, average mw
    ---------------------------------------------------------------------------------------------"""
    n_mw = 0
    mw_sum = 0.0
    for mw in mwlist:
        n_mw += 1
        mw_sum += mw

    average = 0
    if n_protein > 0:
        average = mw_sum / n_mw
    else:
        sys.stderr.write('Fasta::mw_average - error calculating mw, no values in list\n')

    return average


def composition_sum(composition):
    """---------------------------------------------------------------------------------------------
    calculates the total composition of each letter.

    :param composition: list of dict of int, keys are sequence characters, values are counts
    :return: dict of int
    ---------------------------------------------------------------------------------------------"""
    n_char = 0
    sum = {}
    for eachprotein in composition:

        for letter in eachprotein:
            n_char += eachprotein[letter]
            try:
                sum[letter] += eachprotein[letter]
            except KeyError:
                sum[letter] = eachprotein[letter]

    return sum


# ==================================================================================================
# testing - do not change the main program
# ==================================================================================================
if __name__ == '__main__':

    filename = 'test.fasta'
    protein = Fasta(filename)

    composition = []
    mweight = []

    n_protein = 0
    n_residue = 0
    while protein.next():

        n_protein += 1
        n_residue += len(protein.sequence)
        sys.stdout.write('{}\t{} residues\n'.format(protein.id, len(protein.sequence)))

        composition.append(protein.composition())
        mweight.append(protein.mw())

    # report
    sys.stdout.write('\n{} proteins read from {}\n'.format(n_protein, filename))
    sys.stdout.write('\taverage molecular weight: {:.2f}\n'.format(mw_average(mweight)))
    sys.stdout.write('\taverage composition (percent):\n')

    sum_comp = composition_sum(composition)
    for aa in sum_comp:
        sys.stdout.write(
            '{:>10}{:>10d}{:>10.3f}\n'.format(aa, sum_comp[aa], 100*sum_comp[aa] / n_residue))

    exit(0)

import sys


# ==================================================================================================
class Fasta:
    """=============================================================================================
    Reads FastA sequences sequentially from a file, see testing section for usage examples
    ============================================================================================="""
    # molecular wt data is a class variable
    aa_mw = { 'A': 100, 'C':200, 'D':300}

    def __init__(self, filename):
        """-----------------------------------------------------------------------------------------
        Fasta object external attributes
            id
            documentation
            sequence
         internal attributes
            filename
            fh - filehandle
            buffer - string with the next line in the file

        :param filename: string, must be a readable file name
        -----------------------------------------------------------------------------------------"""
        self.filename = filename
        self.fh = None
        self.id = ''
        self.doc = ''
        self.seq = ''
        self.buffer = ''

        try:
            self.fh = open(filename, 'r')
        except IOError:
            sys.stderr.write('Fasta::__init__ - error opening file ({})\n'.format(filename))

        # the file may have blank lines so read down to a non-blank line beginning with >
        while not self.buffer.startswith('>'):
            self.buffer = self.fh.readline().rstrip()

    def next(self):
        """-----------------------------------------------------------------------------------------
        Read and parse the next entry in the file.  Return False if no entry is read (e.g., end of
        file has been reached)

        :return: logical
        -----------------------------------------------------------------------------------------"""

        self.id = ''
        self.documentation = ''
        self.sequence = ''

        # parse the title line which is held in buffer
        id = ''
        doc = ''
        try:
            id, doc = self.buffer.split(' ', maxsplit=1)
        except ValueError:
            # if there is no documentation
            id = self.buffer.rstrip()

        self.id = id.lstrip('>')
        self.documentation = doc

        # read the sequence and concatenate until a >
        while self.buffer:
            self.buffer = self.fh.readline()
            if self.buffer.startswith('>'):
                break

            self.sequence += self.buffer.rstrip(' \n')

        if self.sequence:
            return True
        else:
            # if sequence is empty, nothing was read == eof
            return False

    def composition(self, force=True):
        """-----------------------------------------------------------------------------------------
        Return a dictionary with the composition of the sequence.  If force is True, the letters are
        forced into uppercase.

        :return: dict, aa residues or bases as keys
        -----------------------------------------------------------------------------------------"""
        count = {}
        for letter in self.sequence:
            if force:
                letter = letter.upper()
            try:
                count[letter] += 1
            except KeyError:
                count[letter] = 1

        return count

    def mw(self):
        """-----------------------------------------------------------------------------------------
        Return the calculated molecular weight.  Unknown amino acids are treated as weight = 0.0.

        :return:
        -----------------------------------------------------------------------------------------"""
        comp = self.composition()
        mw = 0
        for aa in comp:
            if aa in Fasta.aa_mw:
                mw += Fasta.aa_mw[aa]

        # add one water for the ends
        return  mw + 18.0


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
# testing
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
    sys.stdout.write('\taverage composition:\n')

    sum_comp = composition_sum(composition)
    for aa in sum_comp:
        sys.stdout.write(
            '{:>10}{:>10d}{:>10.3f}\n'.format(aa, sum_comp[aa], 100*sum_comp[aa] / n_residue))

    exit(0)

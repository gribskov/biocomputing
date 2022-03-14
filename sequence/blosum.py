"""=================================================================================================
perform blosum calculation on example data
could be extended to full calculation

Michael Gribskov     09 March 2022
================================================================================================="""
from math import log


def composition_ref(select, frequency=True):
    """---------------------------------------------------------------------------------------------
    Return composition of protein database

    Data from:http://Martz.MolviZ.Org,  April, 2020 by Eric Martz,

    Sources

    "creighton": Thomas E. Creighton (EMBO Heidelberg), Proteins, Structures and Molecular Properties, 2nd ed. 1993,
    W. H. Freeman and Co. Page 5:
    Table of amino acid residue masses and frequencies in 1,021 unrelated proteins of known sequence

    "trinquier": Trinquier G, Sanejouand YH. Which effective property of amino acids is best preserved by the genetic
    code?. Protein Eng. 1998;11(3):153–169. doi:10.1093/protein/11.3.153
    https://pubmed.ncbi.nlm.nih.gov/9613840/
    values from this table: http://proteinsandproteomics.org/content/free/tables_1/table08.pdf
    "105,990 sequences in the nonredundant OWL protein database (release 26.0 e) "

    "carugo": Carugo, Oliviero. “Amino acid composition and protein dimension.”
    Protein science 17,12 (2008): 2187-91. doi:10.1110/ps.037762.108
    Values read from Figure 2 for lengths of 200 residues.
    Derived from uniref50 FASTA file from UniProt, sequence identity <50%. About 550,000 sequences, which give several
    thousand sequences at each length.
    mrg: i assume carugo100 is for proteins > 100 residues

    "swiss2022":UniProtKB/Swiss-Prot protein knowledgebase release 2022_01 statistics

    :param select: string, key of composition to select
    :param frequency: boolean, if True convert to frequency
    :return: dict, keys are amino acids (one letter code)
    ---------------------------------------------------------------------------------------------"""
    data = {'creighton': {'A': 8.3, 'C': 1.7, 'D': 5.3, 'E': 6.2, 'F': 3.9,
                          'G': 7.2, 'H': 2.2, 'I': 5.2, 'K': 5.7, 'L': 9.0,
                          'M': 2.4, 'N': 4.4, 'P': 5.1, 'Q': 4.0, 'R': 5.7,
                          'S': 6.9, 'T': 5.8, 'V': 6.6, 'W': 1.3, 'Y': 3.2},
            'trinquier': {'A': 7.5, 'C': 1.8, 'D': 5.2, 'E': 6.3, 'F': 3.9,
                          'G': 7.1, 'H': 2.2, 'I': 5.5, 'K': 5.8, 'L': 9.1,
                          'M': 2.8, 'N': 4.6, 'P': 5.1, 'Q': 4.1, 'R': 5.2,
                          'S': 7.4, 'T': 6.0, 'V': 6.5, 'W': 1.3, 'Y': 3.3},
            'carugo200': {'A': 8.7, 'C': 1.5, 'D': 5.1, 'E': 6.2, 'F': 4.0,
                          'G': 6.8, 'H': 2.2, 'I': 5.7, 'K': 5.3, 'L': 9.8,
                          'M': 2.4, 'N': 3.9, 'P': 5.0, 'Q': 3.9, 'R': 6.2,
                          'S': 7.0, 'T': 5.3, 'V': 6.5, 'W': 1.3, 'Y': 2.9},
            'carugo100': {'A': 8.5, 'C': 1.9, 'D': 4.7, 'E': 6.0, 'F': 4.0,
                          'G': 6.5, 'H': 2.3, 'I': 5.8, 'K': 5.6, 'L': 9.6,
                          'M': 2.8, 'N': 3.7, 'P': 4.9, 'Q': 3.9, 'R': 6.7,
                          'S': 7.1, 'T': 5.2, 'V': 6.4, 'W': 1.4, 'Y': 2.7},
            'swiss2022': {'A': 8.25, 'C': 1.38, 'D': 5.46, 'E': 6.72, 'F': 3.86,
                          'G': 7.07, 'H': 2.27, 'I': 5.91, 'K': 5.80, 'L': 9.65,
                          'M': 2.41, 'N': 4.06, 'P': 4.74, 'Q': 3.93, 'R': 5.53,
                          'S': 6.64, 'T': 5.35, 'V': 6.86, 'W': 1.10, 'Y': 2.92, }
            }

    sum = 0
    if frequency:
        for aa in data[select]:
            sum += data[select][aa]

        f = {}
        for aa in data[select]:
            f[aa] = data[select][aa] / sum

        return f

    else:
        return data[select]

class transition



def identity(block):
    """---------------------------------------------------------------------------------------------
    Return identity vector of sequence weights

    :param block: dict, key=name of sequence, value = sequence
    :return: dict, key=name of sequence, value = 1
    ---------------------------------------------------------------------------------------------"""
    # weight = {}
    # for seq in block:
    #     weight[seq] = 1

    return {k: 1 for k in block.keys()}


def block_length(block):
    """---------------------------------------------------------------------------------------------
    Return length of sequence in block. Length is the minimum sequence length of all the sequences

    :param block: dict, key=name of sequence, value = sequence
    :return: int, block length
    ---------------------------------------------------------------------------------------------"""
    length = 10000000
    for seq in block:
        length = min(length, len(block[seq]))

    return length


def composition_count(block, blocklen, weight):
    """---------------------------------------------------------------------------------------------
    Calculate the composition of a block

    :param block: dict, key=name of sequence, value = sequence
    :return: 2 dicts, first is composition by column, second is overall composition
    ---------------------------------------------------------------------------------------------"""
    comp = []
    overall = {}
    for pos in range(blocklen):
        comp.append({})
        for seq in block:
            try:
                comp[pos][block[seq][pos]] += weight[seq]
            except KeyError:
                comp[pos][block[seq][pos]] = weight[seq]

    for c in comp:
        for letter in c.keys():
            try:
                overall[letter] += c[letter]
            except KeyError:
                overall[letter] = c[letter]

    return comp, overall


def transition_count(composition):
    """---------------------------------------------------------------------------------------------
    calculate transition counts across all columns of composition

    :param composition: list of dict, count of letters at each position
    :return: dict of dict, symmetric matrix of letter transitions (float)
    ---------------------------------------------------------------------------------------------"""
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    transition = {letter: {letter: 0 for letter in alphabet} for letter in alphabet}

    for c in composition:
        letters = sorted(list(c.keys()))
        for i in range(len(letters)):
            for j in range(i, len(letters)):
                ci = letters[i]
                cj = letters[j]
                if ci == cj:
                    transition[ci][ci] += (c[ci] * (c[ci] - 1)) / 2
                else:
                    transition[ci][cj] += c[ci] * c[cj]
                    # transition[cj][ci] += c[ci] * c[cj]

    return transition


def composition_addprior(prior, composition):
    """---------------------------------------------------------------------------------------------
    add the frequencies/counts in prior to each position of the composition

    :param prior: dict, keys are amino acids
    :param composition: list of dict, positions dictionary of composition

    :return: count at last position (all position should be the same)
    ---------------------------------------------------------------------------------------------"""
    for c in composition:
        sum = 0
        for letter in c:
            c[letter] += prior[letter]
            sum += c[letter]

    return sum


def transition_background(expected, n):
    """---------------------------------------------------------------------------------------------
    calculate the expected transition frequences based on the expected count (background).
    t(i,i) = exp(i)**2
    t(i,j) = s * exp(i) * exp(j)

    :param expected: dict, expected frequences indexed by aa
    :parameter n: int, expected number of transitions
    :return: dict of dict, expected transition frequencies between letters
    ---------------------------------------------------------------------------------------------"""
    alphabet = ''.join(expected.keys())
    transition = {letter: {letter: 0 for letter in alphabet} for letter in alphabet}

    for i in range(len(alphabet)):
        aa_i = alphabet[i]
        for j in range(i, len(alphabet)):
            aa_j = alphabet[j]
            if i == j:
                transition[aa_i][aa_j] = n * expected[aa_i] ** 2
            else:
                transition[aa_i][aa_j] = n * 2 * expected[aa_i] * expected[aa_j]

    return transition


def tabular(table, format):
    """---------------------------------------------------------------------------------------------

    values are printed with f format, labels with s format

    :param table: dict of dict, keys are aa, values assumed to be float
    :param float: string, format for values, e.g., 5.2
    :return:
    ---------------------------------------------------------------------------------------------"""
    out = ''

    if format.find('.'):
        width, precision = format.split('.')
        width = int(width)
    fs = f'{width}s'
    fv = format + 'f'

    # column labels
    column = f'{" ":>{fs}}'
    for aa in table.keys():
        column += f'{aa:>{fs}}'

    out += f'{column}\n'

    for i in table:
        out += f'{i:<{fs}}'
        for j in table:
            try:
                out += f'{table[i][j]:{fv}}'
            except TypeError:
                # capture the Nones
                a = 'NA'
                out += f'{a:>{fs}}'

        out += '\n'

    return out


def logodds(foreground, background, base=2):
    """---------------------------------------------------------------------------------------------
    calculate the symmetric log odds values

    :param foreground: dict of dict, keys=aa
    :param background: dict of dict, keys=aa
    :return: dict of dict
    ---------------------------------------------------------------------------------------------"""
    alpha = ''.join(foreground.keys())
    logodds = {a1: {a2: None for a2 in alpha} for a1 in alpha}

    for i in range(len(alpha)):
        a1 = alpha[i]
        for j in range(i, len(alpha)):
            a2 = alpha[j]
            try:
                odds = foreground[a1][a2] / background[a1][a2]
                if odds == 0:
                    continue
                else:
                    # make symmetric
                    logodds[a1][a2] = log(odds, base)
                    logodds[a2][a1] = log(odds, base)

            except ZeroDivisionError:
                continue

    return logodds


def transition_subset(trans, subset):
    """---------------------------------------------------------------------------------------------
    return a subset of the transition matrix with only the rows and columns in the string subset
    
    :param trans: dict of dict, transition values
    :param subset: string, list of desired keys (amino acids)
    :return: dict of dict, smaller transition matrix
    ---------------------------------------------------------------------------------------------"""
    smaller = {a1: {a2: 0 for a2 in subset} for a1 in subset}

    for a1 in subset:
        for a2 in subset:
            smaller[a1][a2] = trans[a1][a2]

    return smaller


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # METHYLTRANSFERASE BI
    block = {'TCMN_STRGA': 'IADLGGGDGWFLAQILRRHPHATGLLMDLPRVA',
             'TCMO_STRGA': 'FVDLGGARGNLAAHLHRAHPHLRATCFDLPEME',
             'ZRP4_MAIZE': 'LVDVGGGIGAAAQAISKAFPHVKCSVLDLAHVV',
             'CHMT_POPTM': 'LVDVGGGTGAVVNTIVSKYPSIKGINFDLPHVI',
             'COMT_EUCGU': 'VVDVGGGTGAVLSMIVAKYPSMKGINFDLPHVI',
             'COMT_MEDSA': 'LVDVGGGTGAVINTIVSKYPTIKGINFDLPHVI',
             'CRTF_RHOSH': 'LMDVGGGTGAFLAAVGRAYPLMELMLFDLPVVA',
             'OMTA_ASPPA': 'VVDVGGGRGHLSRRVSQKHPHLRFIVQDLPAVI', }

    nseq = len(block)
    weight = identity(block)
    blocklen = block_length(block)

    exp = composition_ref('swiss2022')
    expected_n = nseq * (nseq - 1) * blocklen / 2
    background = transition_background(exp, expected_n)

    sum = 0
    for a1 in 'ACDEFGHIKLMNPQRSTVWY':
        for a2 in 'ACDEFGHIKLMNPQRSTVWY':
            sum += background[a1][a2]
    print(f'background total: {sum}')

    composition, count = composition_count(block, blocklen, weight)
    composition_addprior(exp, composition)
    trans = transition_count(composition)

    lo = logodds(trans, background)

    print(tabular(trans, '5.1'))
    print(tabular(background, '5.1'))
    print(tabular(lo, '5.1'))

    # strans = transition_subset(trans,'FILV')
    # print(tabular(strans, '5.1'))
    # sbackground = transition_subset(background, 'FILV')
    # print(tabular(sbackground, '5.1'))
    # slogodds = transition_subset(logodds, 'FILV')
    # print(tabular(slogodds, '5.1'))

    # for column 1
    blocklen = 1
    expected_n = nseq * (nseq - 1) * blocklen / 2
    exp_sub = {}
    for aa in 'FILV':
        exp_sub[aa] = exp[aa]

    background = transition_background(exp_sub, expected_n)
    composition, count = composition_count(block, blocklen, weight)
    trans = transition_count(composition)
    for a1 in background.keys():
        for a2 in background.keys():
            background[a1][a2] *= expected_n


    strans = transition_subset(trans, 'FILV')
    sbackground = transition_subset(background, 'FILV')
    lo = logodds(strans, sbackground)
    slo = transition_subset(lo, 'FILV')

    print('\ntransition')
    print(tabular(strans, '7.3'))
    print('\nbackground')
    print(tabular(sbackground, '7.3'))
    print('\nlog-odds')
    print(tabular(slo, '7.3'))

    exit(0)

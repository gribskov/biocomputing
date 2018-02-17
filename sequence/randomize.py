"""=================================================================================================
randomize.py

functions for sequence randomization.  Functions return a string with the desired characteristics.

Synopsis
    import randomize

Michael Gribskov     17 February 2018
================================================================================================="""


def composition(seq='', k=1):
    """---------------------------------------------------------------------------------------------
    Count the kmer words in sequnce and return the number of words and a dictionary with the count
    of each

    :param seq: string
    :param k: word length
    :return: n_words, composition
    ---------------------------------------------------------------------------------------------"""
    n = len(seq) - k + 1
    comp = {}
    for pos in range(n):
        try:
            comp[seq[pos:pos + k]] += 1
        except KeyError:
            comp[seq[pos:pos + k]] = 1

    return n, comp


def frequency(comp, n=0):
    """---------------------------------------------------------------------------------------------
    Calculate the frequency of each word given the total number of words and observed counts.
    :param comp:
    :return:
    ---------------------------------------------------------------------------------------------"""
    if n == 0:
        # number of words not provided, count the words
        for word in comp:
            n += comp[word]
    if n <= 0:
        print('\nrandomize::frequency - number of words <= 0')
        exit(1)

    freq = {}
    for word in comp:
        freq[word] = comp[word] / n

    return freq


def withReplace(seq='', k='1'):
    """---------------------------------------------------------------------------------------------
    Sample with replacement.  Composition is not guaranteed to be the same as the original
    sequence but should be close.

    :param sequence: string, sequence to randomize
    :param k: word size, default=1
    :return: string
    ---------------------------------------------------------------------------------------------"""
    n, comp = composition(seq, k)


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    seq1 = 'AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'

    for k in range(1, 4):
        print('\ncomposition - {}mer'.format(k))
        print('    sequence:', seq1)
        nword, comp = composition(seq1, k)
        print('    n words:', nword)
        print('    composition:', comp)

        freq = frequency({}, 0)
        print('frequency:', freq)

    exit(0)

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

    exit(0)

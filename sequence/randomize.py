"""=================================================================================================
randomize.py

functions for sequence randomization.  Functions return a string with the desired characteristics.

Synopsis
    import randomize

Michael Gribskov     17 February 2018
================================================================================================="""
import random


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
    Calculate the frequency (probability) of each word given the total number of words and observed
    counts.

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


def tabular(*tables):
    """---------------------------------------------------------------------------------------------
    print tables in horizonatal tabular format
    TODO: aggregate keys across all tables
    TODO: specify indent
    TODO: automatically determine field width

    :param tables: 1 or more composition/freuency tables
    :return: None
    ---------------------------------------------------------------------------------------------"""
    indent = 4
    field = 7
    precision = 2

    space = ' ' * indent
    format_s = '{{:>{}}}'.format(field)
    format_f = '{{:{}.{}f}}'.format(field, precision)

    for word in sorted(tables[0]):
        print(format_s.format(word), end='')
    print()

    for table in tables:
        for word in sorted(table):
            print('{:7.2f}'.format(table[word]), end='')
        print()

    return None


def withReplace(seq='', k=1):
    """---------------------------------------------------------------------------------------------
    Sample with replacement.  Composition is not guaranteed to be the same as the original
    sequence but should be close.

    TODO: only works for k=1

    :param sequence: string, sequence to randomize
    :param k: word size, default=1
    :return: string
    ---------------------------------------------------------------------------------------------"""
    n, comp = composition(seq, k)
    freq = frequency(comp)

    # generate cumulative frequencies (thresholds)
    threshold = {}
    sum = 0
    for word in sorted(freq):
        sum += freq[word]
        threshold[word] = sum

    rseq = ''
    found = ''
    for i in range(len(seq)):
        r = random.random()
        for word in sorted(freq):
            if r <= threshold[word]:
                found = word
                break
            found = word
        rseq += found

    return rseq


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    seq1 = 'AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'

    for k in range(1, 4):
        print('\ncomposition - {}mer'.format(k))
        nword, comp = composition(seq1, k)
        freq = frequency(comp, 0)

        print()
        tabular(comp, freq)

    seq2 = 'A' * 100 + 'C' * 100 + 'G' * 100 + 'T' * 100
    n, comp = composition(seq2)
    freq = frequency(comp)

    rseq = withReplace(seq2)
    nr, rcomp = composition(rseq)
    rfreq = frequency(rcomp)
    print()
    tabular(rfreq, freq)

    exit(0)

"""=================================================================================================
wordmatch

word related matching for sequences

identities in window
total score withing window using scoring table

TODO: conversion between diaganl, offset and sequence coordinates


Michael Gribskov     17 June 2020
================================================================================================="""
from sequence.fasta import Fasta


class Match():
    """=============================================================================================

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.diagonal = []

    def identity(self, s1, s2):
        """-----------------------------------------------------------------------------------------

        :param s1: Fasta, sequence 1
        :param s2: Fasta, sequence 2
        :return: int, number of matches
        -----------------------------------------------------------------------------------------"""
        l1 = len(s1.seq)
        l2 = len(s2.seq)

        nmatch = 0
        for diag in range(l1 + l2 - 1):
            pos1 = max(diag - l2 + 1, 0)
            pos2 = max(l2 - diag - 1, 0)
            n = min(l1 - pos1, l2 - pos2)
            self.diagonal.append([])
            for offset in range(n):
                # print('s1:{}     s2:{}     diag: {}    n: {}'.format(pos1, pos2, diag, n))
                if s1.seq[pos1] == s2.seq[pos2]:
                    self.diagonal[diag].append(offset)
                    nmatch += 1
                pos1 += 1
                pos2 += 1

        return nmatch

    # --------------------------------------------------------------------------------------------------


# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    print('\ntest 0: identity matching')
    print('\texpect 7 matches\n')
    fasta = Fasta()
    fasta.id = 'test0'
    fasta.doc = '5 letter DNA test'
    fasta.seq = 'ACAGT'
    print('{}\n'.format(fasta.format()))

    match = Match()
    nmatch = match.identity(fasta, fasta)
    print('matches: {}'.format(nmatch))

    print('\ntest 1: identity matching, unequal length sequences')
    print('\texpect 11 matches\n')
    fasta1 = Fasta()
    fasta1.id = 'test1.1'
    fasta1.doc = '5 letter DNA test'
    fasta1.seq = 'ACAGT'
    print('{}\n'.format(fasta1.format()))

    fasta2 = Fasta()
    fasta2.id = 'test1.2'
    fasta2.doc = '7 letter DNA test'
    fasta2.seq = 'ACAGTAA'
    print('{}\n'.format(fasta2.format()))

    match = Match()
    nmatch = match.identity(fasta1, fasta2)
    print('matches: {}'.format(nmatch))

    exit(0)

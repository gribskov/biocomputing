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
        self.s1 = None
        self.s2 = None

    def identityPos(self):
        """-----------------------------------------------------------------------------------------
        MAach identical characters and return as a list of diagonals with matching positions

        :return: int, number of matches
        -----------------------------------------------------------------------------------------"""
        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)

        nmatch = 0
        for diag in range(l1 + l2 - 1):
            pos1 = max(diag - l2 + 1, 0)
            pos2 = max(l2 - diag - 1, 0)
            n = min(l1 - pos1, l2 - pos2)
            self.diagonal.append([])
            for offset in range(n):
                # print('s1:{}     s2:{}     diag: {}    n: {}'.format(pos1, pos2, diag, n))
                if s1[pos1] == s2[pos2]:
                    self.diagonal[diag].append(offset)
                    nmatch += 1
                pos1 += 1
                pos2 += 1

        return nmatch

    def identity(self):
        """-----------------------------------------------------------------------------------------
        Identify matching regions on diagonals with run length encoding.

        :return: int, number of matches
        -----------------------------------------------------------------------------------------"""
        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)

        nmatch = 0
        for diag in range(l1 + l2 - 1):
            pos1 = max(diag - l2 + 1, 0)
            pos2 = max(l2 - diag - 1, 0)
            n = min(l1 - pos1, l2 - pos2)
            self.diagonal.append([])
            run = 0
            for offset in range(n):
                # print('s1:{}     s2:{}     diag: {}    n: {}'.format(pos1, pos2, diag, n))
                if s1[pos1] == s2[pos2]:
                    run += 1
                    nmatch += 1
                elif run:
                    self.diagonal[diag].append([offset, run])
                    run = 0

                pos1 += 1
                pos2 += 1

            if run:
                self.diagonal[diag].append([offset, run])

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
    match.s1 = fasta
    match.s2 = fasta
    nmatch = match.identityPos()
    print('matches: {}'.format(nmatch))

    print('\ntest 1: identity matching, unequal length sequences')
    print('\texpect 11 matches\n')
    match = Match()

    fasta1 = Fasta()
    fasta1.id = 'test1.1'
    fasta1.doc = '5 letter DNA test'
    fasta1.seq = 'ACAGT'
    match.s1 = fasta1
    print('{}\n'.format(match.s1.format()))

    fasta2 = Fasta()
    fasta2.id = 'test1.2'
    fasta2.doc = '7 letter DNA test'
    fasta2.seq = 'ACAGTAA'
    match.s2 = fasta2
    print('{}\n'.format(match.s2.format()))

    nmatch = match.identity()
    print('matches: {}'.format(nmatch))

    exit(0)

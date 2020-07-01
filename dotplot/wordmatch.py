"""=================================================================================================
wordmatch

word related matching for sequences based on exact matching of letters.

identities in a window of defined length
filters by run length
filters by presence in window of minimum score

use dotplot.plotter to display

Michael Gribskov     17 June 2020
================================================================================================="""
from sequence.fasta import Fasta


class Base():
    """=============================================================================================
    Base class for matching codes.  Contains methods useful for all classes

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Constructor for Base object.  All match classes require sequences.
        -----------------------------------------------------------------------------------------"""
        self.s1 = None
        self.s2 = None


    def rle2coord(self):
        """-----------------------------------------------------------------------------------------
        Return a list of beginning and ending positions of each run.  List is a list of four
        coordinates for each run [s1begin, s1end, s2begin, s2end]

        :return:
        -----------------------------------------------------------------------------------------"""
        coord = []

        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)

        for diag in range(len(self.diagonal)):

            for offset, length in self.diagonal[diag]:
                end1 = max(diag - l2 + 1, 0) + offset
                end2 = max(l2 - diag - 1, 0) + offset
                beg1 = end1 - length + 1
                beg2 = end2 - length + 1
                coord.append([beg1, end1, beg2, end2])

        return coord

    @staticmethod
    def diagLenBegin(diag, l1, l2):
        """-----------------------------------------------------------------------------------------
        Calculates the length of diagonal diag and the beginning position of the diagonal in
        each sequence

        :param diag: int, diagonal number
        :param l1: int, length of sequence 1
        :param l2: int, lenght of sequence 2
        :return: int (diagonal length), int (seq1 begin), int (seq2 begin)
        -----------------------------------------------------------------------------------------"""
        pos1 = max(diag - l2 + 1, 0)
        pos2 = max(l2 - diag - 1, 0)
        diaglen = min(l1 - pos1, l2 - pos2)

        return diaglen, pos1, pos2


class Match(Base):
    """=============================================================================================
    A container for the sequences and the match between them.  The match is stored in diagonal as
    run length encoded pairs of [offser, length] where offset is the position of the end of the
    run along the diagona (zero origin)

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Constructor for Match object.  Match object is an RLE list of diagonals with matching areas.
        each run of matches is encoded as [offset, length] in self.diagonal

        -----------------------------------------------------------------------------------------"""
        self.diagonal = []

    def identityPos(self):
        """-----------------------------------------------------------------------------------------
        Match identical characters and return as a list of diagonals with matching positions

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
            diaglen, pos1, pos2 = Match.diagLenBegin(diag, l1, l2)
            # pos1 = max(diag - l2 + 1, 0)
            # pos2 = max(l2 - diag - 1, 0)
            # n = min(l1 - pos1, l2 - pos2)
            self.diagonal.append([])
            runlen = 0
            for offset in range(diaglen):
                # print('s1:{}     s2:{}     diag: {}    n: {}'.format(pos1, pos2, diag, n))
                if s1[pos1] == s2[pos2]:
                    runlen += 1
                    nmatch += 1
                elif runlen:
                    # no match, end of run. this offset is one past the end of the run
                    self.diagonal[diag].append([offset - 1, runlen])
                    runlen = 0

                pos1 += 1
                pos2 += 1

            if runlen:
                # end of a run at end of diagonal (do not subtract because the end position is
                # the true end)
                self.diagonal[diag].append([offset, runlen])

        return nmatch

    def filterByLen(self, n):
        """-----------------------------------------------------------------------------------------
        remove runs shorter than n

        :param n: int, minimum run length
        :return: int, number of runs
        -----------------------------------------------------------------------------------------"""
        nrun = 0
        for d in self.diagonal:

            l = len(d)
            while l > 0:
                # work from end to beginning because deleting elements changes the length of the
                # list
                l -= 1
                if d[l][1] < n:
                    del (d[l])
                else:
                    nrun += 1

        return nrun

    def filterByWindowCount(self, window, count):
        """-----------------------------------------------------------------------------------------
        For a window of length w, keep the match if the window contains at least count matches.
        This is different that the conventional practice of just marking a single position at the
        center of the window.

        :param window: int, window length
        :param count: int, minimum count to keep
        :return: int, number of runs
        -----------------------------------------------------------------------------------------"""
        l1 = len(self.s1.seq)
        l2 = len(self.s2.seq)

        diagonal = self.diagonal
        for d in range(len(diagonal)):
            diaglen, begin1, begin2 = Match.diagLenBegin(d, l1, l2)

            if diaglen < window:
                # skip  and filter diagonals shorter than window length
                diagonal[d] = []
                continue
            if not len(diagonal[d]):
                # skip diagonals with no runs at all
                continue

            tmpdiag = [0 for i in range(diaglen)]
            for rend, rlen in diagonal[d]:

                # copy runs into temporary diagonal
                for pos in range(rend - rlen + 1, rend + 1):
                    tmpdiag[pos] = 1

            # count the first window
            sum = 0
            for pos in range(window):
                sum += tmpdiag[pos]

            wend = 0
            if sum >= count:
                wend = window

            for pos in range(diaglen - window):

                if sum >= count:
                    wend = pos + window

                olddiag = tmpdiag[pos]
                if pos < wend:
                    tmpdiag[pos] += 1

                sum += tmpdiag[pos + window]
                sum -= olddiag

            # add one to positions in the last positive window
            if sum >= count:
                wend = pos + window + 1

            for pos in range(diaglen - window, diaglen):
                if pos < wend:
                    tmpdiag[pos] += 1
                else:
                    break

            # now find the runs, all matches in a sufficiently high scoring window have been
            # incremented so test for > 1 instead of >=

            runlen = 0
            nmatch = 0
            filtered = []
            for pos in range(diaglen):
                # print('s1:{}     s2:{}     diag: {}    n: {}'.format(pos1, pos2, diag, n))
                if tmpdiag[pos] > 1:
                    runlen += 1
                    nmatch += 1
                elif runlen:
                    # no match, end of run. this offset is one past the end of the run
                    filtered.append([pos - 1, runlen])
                    runlen = 0

            if runlen:
                # end of a run at end of diagonal (do not subtract because the end position is
                # the true end)
                filtered.append([pos, runlen])

            diagonal[d] = filtered

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
    coord = match.rle2coord()

    exit(0)

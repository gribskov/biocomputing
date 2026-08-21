"""=================================================================================================
just a trial of using an object to represent the cell and pointers.  since it is a local alignment
only non'zero cells need to be represented.  maybe this will make it sparse enough to be worthwhile

Michael Gribskov     04 August 2020
================================================================================================="""
import sys
from math import log10
import random
from scipy import stats
import matplotlib.pyplot as plt
from sequence.fasta import Fasta
from sequence.score import Score


class Cell:
    """=============================================================================================

    ============================================================================================="""
    count = 0

    def __init__(self, x=None, y=None, p=None):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.n = Cell.count
        Cell.count += 1
        self.score = 0
        self.p = []
        self.xy = []
        if x and y:
            self.xy = [x, y]
        if p:
            self.p = p


class Alignment(Score):
    """---------------------------------------------------------------------------------------------

    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        Score.__init__(self)
        self.s1 = None
        self.s2 = None
        self.i1 = None
        self.i2 = None
        self.score = None

    def seqToInt(self):
        """-----------------------------------------------------------------------------------------
        Convert sequence strings to an integer arrays and stores in object.  An integer array is
        more convenient for direct lookups in the scoring table than a string

        :return: int, int length of sequence lists
        -----------------------------------------------------------------------------------------"""
        a2i = self.a2i

        self.i1 = [a2i[c] for c in self.s1]
        self.i2 = [a2i[c] for c in self.s2]

        return len(self.i1), len(self.i2)

    def globalBrute(self, open, extend, nogap=False):
        """-----------------------------------------------------------------------------------------
        Local alignment score only.
        s1 is the horizontal sequence and s2 is the vertical sequence.  this makes s2 the row
        index and s1 the column index.

        :param open: float, gap opening penalty
        :param extend: float, gap extension penalty
        :return:
        -----------------------------------------------------------------------------------------"""
        cmp = self.table
        i1 = self.i1
        i2 = self.i2
        l1 = len(i1)
        l2 = len(i2)

        edge = Cell()  # a dummy cell for edges
        edge.score = self.min + open + l1 * extend

        score = [[Cell(i, j) for i in i1] for j in i2]
        self.score = score
        bestrow = Cell()
        bestcol = [Cell() for i in i1]

        score[0][0].score = cmp[i2[0]][i1[0]]
        gap = open
        for ipos in range(l1):
            bestcol[ipos].score = edge.score
            bestcol[ipos].p = []
            gap += extend

        diag = Cell()

        jpos = 0
        vgap = 0
        diag.score = 0
        for j in i2:
            bestrow.p = []
            bestrow.score = edge.score

            ipos = 0
            for i in i1:
                previous = max(diag.score, bestcol[ipos].score, bestrow.score)
                cell = cmp[j][i] + previous

                for dir in (diag, bestrow, bestcol[ipos]):
                    if dir.score == previous:
                        # set pointers for all directions
                        score[jpos][ipos].p.append(dir)

                # update best row and column values
                if diag.score + open > bestrow.score + extend:
                    # what if scores are equal? some paths missed
                    bestrow.p = diag.p
                    bestrow.score = diag.score + open
                else:
                    bestrow.score += extend

                if diag.score + open > bestcol[ipos].score + extend:
                    bestcol[ipos].p = diag.p
                    bestcol[ipos].score = diag.score + open
                else:
                    bestcol[ipos].score += extend

                # diagonal score for next cell
                if jpos > 0:
                    diag.p = score[jpos - 1][ipos]
                    diag.score = score[jpos - 1][ipos].score
                else:
                    diag.score = edge.score

                score[jpos][ipos].score = cell
                ipos += 1

            # end of loop over columns

            # special case for first cell in each row, best previous is always a column gap
            if jpos:
                vgap += extend
            else:
                vgap = open
            bestcol[0].score = vgap

            diag.score = edge.score
            diag.p = []

            jpos += 1

            # end of loop over rows

        # add the end gap penalties
        if not nogap:
            gap = open
            jpos = len(i2) - 1
            for ipos in range(len(i1) - 2, -1, -1):
                score[jpos][ipos].score += gap
                gap += extend

            gap = open
            ipos = len(i1) - 1
            for jpos in range(len(i2) - 2, -1, -1):
                score[jpos][ipos].score += gap
                gap += extend

        scoremax = 1
        posmax = [1, 1]
        return scoremax, posmax

    def localBrute(self, open, extend):
        """-----------------------------------------------------------------------------------------
        Local alignment score only.
        s1 is the horizontal sequence and s2 is the vertical sequence.  this makes s2 the row
        index and s1 the column index.

        :param open: float, gap opening penalty
        :param extend: float, gap extension penalty
        :return:
        -----------------------------------------------------------------------------------------"""
        cmp = self.table
        i1 = self.i1
        i2 = self.i2
        l1 = len(i1)
        l2 = len(i2)

        edge = Cell()  # a dummy cell for the edges
        edge.score = 0

        # set up scoring matrix size l2 * l1, and create x,y position labels
        score = [Cell() for i in range(l2 * l1)]
        self.score = score
        # first row, previous is edge
        for i in range(l1):
            score[i].xy = [i, 0]
            # score[i].p = [edge]

        x = y = 0
        for i in range(l1, len(score)):
            c = score[i]
            if i % l1:
                x += 1
            else:
                # first cell in row
                y += 1
                x = 0
                c.p = [edge]
                c.score = max(0, cmp[i1[x]][i2[y]])
            c.xy = [x, y]

        scoremax = 0
        posmax = []

        # auxiliary storage: 1 pointer for the best gapped value in the previous row (y-1, x:0..-1)
        # 1 pointer for the best gapped value in each column. use a cell object for the pointers
        xgap = Cell()
        ygap = [Cell() for i in range(l1)]

        # first row
        x = y = 0
        for c in self.score[0:l1]:
            c.score = max(0, cmp[i1[x]][i2[y]])
            if c.score >= scoremax:
                scoremax = c.score
                if c.score == scoremax:
                    posmax += [c]
                else:
                    posmax = [c]
            # if c.score + open > 0:
            # ygap[x].p = [edge]
            ygap[x].score = 0
            x += 1

        for c in self.score[l1:]:
            x = c.xy[0]
            y = c.xy[1]
            if x == 0:
                # left edge cell, diag, xgap, and ygap[x-1] undefined
                bestscore = max(0, cmp[i1[0]][i2[y]])
                xgap.score = 0
                # xgap.p = [edge]

            else:
                # internal cell
                diag = score[c.n - l1 - 2]
                bestprevscore = max(diag.score, xgap.score, ygap[x-1].score)
                if diag.score == bestprevscore:
                    c.p += [diag]
                if xgap.score == bestprevscore:
                    c.p += xgap.p
                if ygap[x-1].score == bestprevscore:
                    c.p += ygap[x-1].p

                c.score = max(0, bestprevscore + cmp[i1[x]][i2[y]])
                if c.score == 0: continue

                if c.score >= scoremax:
                    if c.score == scoremax:
                        posmax += [c]
                    else:
                        posmax = [c]
                    scoremax = c.score

                # update gap pointers
                testdiag = diag.score + open
                testx = xgap.score + extend
                if  testdiag > testx:
                    xgap.score = testdiag
                    xgap.p = [diag]
                elif testdiag == testx:
                    xgap.score = testx
                    xgap.p += [diag]
                else:
                    xgap.score = testx

                testy =  ygap[x-1].score + extend
                if testdiag > testy:
                    ygap[x - 1].score = testdiag
                    ygap[x - 1].p = [diag]
                elif testdiag == testy:
                    ygap[x - 1].score = testy
                    ygap[x - 1].p += [diag]
                else:
                    ygap[x-1].score = testy

                x += 1


        return scoremax, posmax

    def trace1(self, pos):
        """-----------------------------------------------------------------------------------------
        Trace back one alignment using the first pointer for each cell

        :param pos: list of 2 int, traceback start position
        :return:
        -----------------------------------------------------------------------------------------"""
        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)
        cmp = self.table
        a2i = self.a2i

        score = self.score
        current = score[pos[0]][pos[1]]

        a1 = ''
        a2 = ''
        rowold = pos[0]
        colold = pos[1]
        while len(current.p) > 0:
            row, col = Alignment.n2pos(l1, current.n)

            for c in range(colold - 1, col, -1):
                a1 += s1[c]
                a2 += '.'

            for r in range(rowold - 1, row, -1):
                a1 += '.'
                a2 += s2[r]

            a1 += s1[col]
            a2 += s2[row]

            if len(current.p):
                current = current.p[0]

            rowold = row
            colold = col

        m = self.matchString(a1, a2)

        return a1[::-1], a2[::-1], m[::-1]

    def traceAllPtr(self, endpts):
        """----------------------------------------------------------------------------------------

        ----------------------------------------------------------------------------------------"""
        print(f'starting traceback')
        s1 = align.s1
        s2 = align.s2
        stack = []
        for c in endpts:
            a1 = s1[c.xy[0]]
            a2 = s2[c.xy[1]]
            stack =[ [c, a1, a2]]

        while stack:
            (c, a1, a2) = stack.pop()
            print(f'pop [{c.xy[0]}, {c.xy[1]}]')
            # cn = val[0]
            # a1 = val[1]
            # a2 = val[2]
            if c.p:

                for cn in c.p:
                    if not cn.xy:
                        print(f'not xy\n{a1[::-1]}\n{a2[::-1]}\n')
                        continue

                    na1 = a1
                    na2 = a2
                    for x in range(c.xy[0] - 1, cn.xy[0], -1):
                        na1 += s1[x]
                        na2 += '.'
                    for y in range(c.xy[1] - 1, cn.xy[1], -1):
                        na1 += '.'
                        na2 += s2[y]

                    na1 += s1[cn.xy[0]]
                    na2 += s2[cn.xy[1]]
                    if cn.xy:
                        stack.append([cn, na1, na2])

            else:
                print(f'not p \n{a1[::-1]}\n{a2[::-1]}\n')

        return

    def traceAll(self, pos):
        """-----------------------------------------------------------------------------------------
        Trace back one alignment using the first pointer for each cell

        :param pos: list of 2 int, traceback start position
        :return:
        -----------------------------------------------------------------------------------------"""
        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)
        cmp = self.table
        a2i = self.a2i

        score = self.score
        stack = []
        a1 = ' ' * (l1 * l2)
        a2 = ' ' * (l1 * l2)
        n = score[pos[0]][pos[1]].n
        nold = score[pos[0]][pos[1]].n
        alen = 0
        stack.append([n, nold, alen])

        save = []
        while stack:
            n, nold, alen = stack.pop()
            row, col = Alignment.n2pos(l1, n)
            rowold, colold = Alignment.n2pos(l1, nold)
            a1 = a1[:alen]
            a2 = a2[:alen]

            for c in range(colold - 1, col, -1):
                a1 += s1[c]
                a2 += '.'
                alen += 1

            for r in range(rowold - 1, row, -1):
                a1 += '.'
                a2 += s2[r]
                alen += 1

            a1 += s1[col]
            a2 += s2[row]
            alen += 1

            # for each path in the pointers of the current cell push on stack
            ptrs = score[row][col].p
            if len(ptrs):
                for p in ptrs:
                    stack.append([p.n, n, alen])

            else:
                # if there are no pointers, it is the end of a path
                save.append([a1[:alen], a2[:alen]])

        return save

    def writeScoreMatrix(self, file, decimal=0, reverse=False, space=2):
        """-----------------------------------------------------------------------------------------
        Write out the score matrix in an aligned table. Could provide scoremax, but then it wouldn't
        work for a sub-table.

        :param file: open filehandle for output, think stdout
        :param decimal:int, number of digits past decimal point
        :return: float, maximum score
        -----------------------------------------------------------------------------------------"""
        # first find the largest value (maximum column width)
        score = self.score
        scoremax = 0.0
        for row in score:
            for col in row:
                scoremax = max(scoremax, col.score)

        fmt = '{{:>{}.{}f}}'.format(len(str(scoremax)) + space, decimal)
        smt = '{{:>{}s}}'.format(len(str(scoremax)) + space)

        if reverse:
            # up, left
            s1 = self.s1.seq
            s2 = self.s2.seq
            i0 = len(s1) - 1
            j0 = len(s2) - 1
            step = -1

            file.write(smt.format(' '))
            i = i0
            while i >= 0:
                file.write(smt.format(s1[i]))
                i += step
            file.write('\n')

            j = j0
            while j >= 0:
                file.write(smt.format(s2[j]))
                i = i0
                while i >= 0:
                    s = score[i][j]
                    file.write(fmt.format(score[j][i].score))
                    i += step
                j += step

                file.write('\n')
        else:
            # down, right
            file.write(smt.format(' '))
            for c in self.s1.seq:
                file.write(smt.format(c))
            file.write('\n')

            s2 = self.s2.seq
            i = 0
            for row in score:
                file.write(smt.format(s2[i]))
                i += 1
                for col in row:
                    file.write(fmt.format(col.score))

                file.write('\n')

        return scoremax

    def matchString(self, a1, a2):
        """-----------------------------------------------------------------------------------------
        creat a string showing the match between the aligned sequences a1 and a2.

        :param a1: string, aligned sequence 1
        :param a2: string, aligned sequence 2
        :return: string, match string
        -----------------------------------------------------------------------------------------"""
        idchar = '|'
        simchar = ':'
        cmp = self.table
        a2i = self.a2i

        match = ''
        for i in range(len(a1)):
            c1 = a1[i]
            c2 = a2[i]
            if c1 == c2:
                match += idchar
            elif c1 == '.' or c2 == '.':
                match += ' '
            elif cmp[a2i[c1]][a2i[c2]] > 0:
                match += simchar
            else:
                match += ' '

        return match

    @staticmethod
    def n2pos(l1, n):
        """-----------------------------------------------------------------------------------------
        return the row and col corresponding to cell n. Assumes the matrix is stored in row major
        order with l2 rows (y) and l1 columns (x)

        :param l1: int, length of sequence 1 (col)
        :param n: int, cell n
        :return: int, int; row, col
        -----------------------------------------------------------------------------------------"""
        return (n - 1) // l1, (n - 1) % l1

    def draw_grid_lines(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        fig, ax = plt.subplots(figsize=(6, 6))
        i_size = len(self.i2)
        j_size = len(self.i1)

        # Draw grid lines
        for x in range(i_size + 1):
            ax.axvline(x, color='lightgray', linestyle='--')
        for y in range(j_size + 1):
            ax.axhline(y, color='lightgray', linestyle='--')

        # Draw connections from cell center to target centers
        for x in range(len(self.i1)):
            cy = -0.75
            cx = x + 0.5
            ax.text(cx, cy + 0.1, self.s1[x], fontsize=10, fontweight='bold', ha='center', va='bottom', color='blue')

        for y in range(len(self.i2)):
            cy = y + 0.25
            cx = -0.5
            ax.text(cx, cy + 0.1, self.s2[y], fontsize=10, fontweight='bold', ha='center', va='bottom', color='blue')


        score = self.score
        for c in score:
            start_x = c.xy[0] + 0.5
            start_y = c.xy[1] + 0.5
            ax.text(start_x, start_y, str(c.score), fontsize=10, fontweight='bold',
                    ha='center', va='center', color='black', zorder=2)
            if c.p:
                for p in c.p:
                    if p.xy:
                        end_x = p.xy[0] + 0.5
                        end_y = p.xy[1] + 0.5
                        # ax.plot([start_x, end_x], [start_y, end_y], color='blue', marker='o')
                        ax.plot([start_x, end_x], [start_y, end_y], color='red', linewidth=0.75)


        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(0, len(self.s1))
        ax.set_ylim(0, len(self.s2))
        cell_label = f"[{self.s2},{self.s1}]"
        # ax.text(cx, cy + 0.1, cell_label, fontsize=8, ha='center', va='bottom', color='blue')
        # ax.set_aspect('equal')
        plt.show()

        return

# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    test_sequence = ['ACTTATCTTAT', 'TATTCTATTCA', 'TGGTATACTAT', 'GATACTATCTA',
                     'AGTATCATATT', 'TTATACTATGG', 'TACTATTTAGAT', 'TTATACTATGA',
                     'TAGATTTATCAT', 'TGGTATACTAT', 'BORROW', 'BORABORA']

    # sequences
    align = Alignment()
    align.s1 = test_sequence[6]
    align.s2 = test_sequence[7]

    # scoring table
    # align.alphabet = 'ACGT'
    # align.identity(pos=3, neg=-3)
    # align.readNCBI('..//tables/alphabet.matrix')
    # align.readNCBI('../../dotplot/table/NUC4.4.matrix')
    align.readNCBI('..//tables/dna4-2.matrix')

    # random.shuffle(align.i1)          # uncomment to test scores for random alignments
    align.seqToInt()
    # bestscore, bestpos = align.globalBrute(-1, -1, nogap=False)
    bestscore, bestpos = align.localBrute(-1, -1)
    align.traceAllPtr(bestpos)
    for c in bestpos:
        print(f'score: {bestscore} at {c.xy}\n')
    # align.writeScoreMatrix(sys.stdout, reverse=False, space=1)
    # alignments = align.traceAll(bestpos)
    align.draw_grid_lines()

    exit(0)

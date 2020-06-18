"""=================================================================================================
dotplot.plotter

Michael Gribskov     18 June 2020
================================================================================================="""
import matplotlib.pyplot as plt


class Plotter():
    """=============================================================================================
    matplotlib plotting class for dotplots
    basically a container for a matplotlib figure and a match object
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.fig = plt.figure()
        self.match = None

    def show(self, *args, **kwargs):
        """-----------------------------------------------------------------------------------------
        Delegate to plt.show().  makes syntax a little easier in application since the object is
        used instead of the class

        :param args:
        :param kwargs:
        :return: True
        -----------------------------------------------------------------------------------------"""
        plt.show(*args, **kwargs)

        return True

    def title(self, titlestr='dotplot'):
        """-----------------------------------------------------------------------------------------
        Sets the title

        :return: True
        -----------------------------------------------------------------------------------------"""
        self.fig.suptitle(titlestr)

        return True

    def setup(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        self.title()
        ax = self.fig.add_subplot(1, 1, 1)

        ax.set_xlim(0, len(self.match.s1.seq) + 1)
        ax.set_ylim(0, len(self.match.s2.seq) + 1)

        ax.set_xlabel('\n'.join([self.match.s1.id, self.match.s1.doc]))
        ax.set_ylabel('\n'.join([self.match.s2.doc, self.match.s2.id]))

        return True

    def lines(self, dither=0.05):
        """-----------------------------------------------------------------------------------------
        plot matches as lines

        :param dither: float, half the length for lines that are only one letter long
        :return:
        -----------------------------------------------------------------------------------------"""
        coord = self.match.rle2coord()
        for line in coord:
            if line[0] == line[1]:
                plt.plot([line[0] + 1 - dither, line[1] + 1 + dither],
                         [line[2] + 1 - dither, line[3] + 1 + dither], color='k')
            else:
                plt.plot([line[0] + 1, line[1] + 1], [line[2] + 1, line[3] + 1], color='k')

        return True

    def dots(self):
        """-----------------------------------------------------------------------------------------
        plot matches as dots

        :return:
        -----------------------------------------------------------------------------------------"""
        coord = self.match.rle2coord()
        for line in coord:
            # print(line)
            if line[0] == line[1]:
                plt.plot(line[0] + 1, line[1] + 1, 'ko')

            y = line[2]
            for x in range(line[0],line[1]+1):
                plt.plot(x + 1, y + 1, 'ko')
                y += 1


        return True
# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    from wordmatch import Match  # for testing only
    from sequence.fasta import Fasta

    print('\ntest 1: identity matching, unequal length sequences')
    print('\texpect 11 matches\n')
    match = Match()

    fasta1 = Fasta()
    fasta1.id = 'test0.1'
    fasta1.doc = '5 letter DNA test'
    fasta1.seq = 'ACAGT'
    match.s1 = fasta1

    fasta2 = Fasta()
    fasta2.id = 'test0.2'
    fasta2.doc = '7 letter DNA test'
    fasta2.seq = 'ACAGTAA'
    match.s2 = fasta2

    nmatch = match.identity()

    plot = Plotter()
    plot.match = match
    plot.setup()
    plot.lines(dither=0.03)
    plot.dots()
    plot.show()

    exit(0)

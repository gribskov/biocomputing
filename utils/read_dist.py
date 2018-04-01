"""=================================================================================================
Plot the distribution of reads on a reference sequence based on mapped reads in a SAM file

SAM format is (all one line, whitespace separated fields)

SRR5295840.120  163     AT1G07250.1     101     44      1S150M  =       347     398     NCT...CAG
#A<...FJF AS:i:300     XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:150        YS:i:300        YT:Z:CP

SRR5295840.120  QNAME   read name
163             FLAG    mapping bit flags
AT1G07250.1     RNAME   reference sequence name
101             POS     leftmost position of mapped read
44              MAPQ    mapping quality
1S150M          CIGAR   alignment
=               RNEXT   name of mate/next read
347             PNEXT   position of mate/next read
398             TLEN    inferred insert size
NCT...CAG       SEQ     sequence
#A<...FJF       QUAL    quality

the remaining columns are application specific
AS:i:300     XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:150        YS:i:300        YT:Z:CP

Michael Gribskov    1 April 2018
================================================================================================="""
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def extend_list(arr, end, init=0):
    """---------------------------------------------------------------------------------------------
    extend the existing list arr by adding indices from the current end of the list to the specified
    end pos (the new last index in list)

    :param arr: list
    :param end: last index to create
    :param init: value to initialize elements with
    :return: int, new list size
    ---------------------------------------------------------------------------------------------"""
    begin = len(arr)
    arr += [init for k in range(begin, end + 1)]

    return len(arr)


def add_cigar(seq, pos, cigar):
    """---------------------------------------------------------------------------------------------
    increment the count of reads at the corresponding positions of sort given the start pos and
    CIGAR string

    :param seq: list corresponding to sequence bases
    :param pos: beginning position of mapped read
    :param cigar: alignment string
    :return: integer, bases incremented (number of M)
    ---------------------------------------------------------------------------------------------"""
    istr = ''
    m = 0
    for char in cigar:
        if char.isdigit():
            istr += char
            continue

        i = int(istr)
        istr = ''
        if char == 'M':
            # matching positions

            # check to make sure seq list is big enough, if not add some more elements
            if pos + i + 1 > len(seq):
                extend_list(seq, pos + i + 10000)

            for j in range(pos, pos + i - 1):
                m += 1
                seq[j] += 1

        elif char not in 'SD':
            # must be a character we don't care about, we'll just ignore these for now
            # I gets dealt with here
            continue

        # only M, S, and D fall through to here
        #  M, S, and D increment the position, S and D do nothing else
        pos += i - 1

        return m

    # end of add_cigar


def condense_seq(seq):
    """---------------------------------------------------------------------------------------------
    condense the seq array by removing consecutive positions with the same value.  This allows
    easier printing. For plotting, each interval ends at the indicated position.

    pos1, val1
    pos2, val2
    ...

    SAM files use a 1 origin, so the first interval is 1 to pos1 with value val1, the second is
    pos1 + 1 to pos2 with val2, etc.

    :param seq: list of int, coverage at each base
    :return: list of tuples, (pos, val)
    ---------------------------------------------------------------------------------------------"""
    compressed = []
    cover = seq[1]
    for pos in range(2, len(seq)):
        if seq[pos] != cover:
            compressed.append((pos - 1, cover))
        cover = seq[pos]

    return compressed


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    map = None
    try:
        map = open(sys.argv[1], 'r')
    except:
        print('unable to open input file ({}'.format(sys.argv[1]))
        exit(1)

    nread = 0
    seq = []

    # assume that the mapped read begins at POS
    # use the CIGAR string to increment counts in the seq array

    bases_mapped = 0
    for line in map:
        if line.startswith('@'):
            # skip header lines
            continue

        field = line.split()
        # print('{}\t{}'.format(field[0], field[8]))
        pos = int(field[3])
        mapq = field[4]
        cigar = field[5]

        # filter some lines
        if mapq == 0 or cigar == '*':
            continue

        nread += 1
        bases_mapped += add_cigar(seq, pos, cigar)

        if nread > 100000:
            break

    print('{} bases mapped'.format(bases_mapped))

    comp_coverage = condense_seq(seq)
    for coverage in comp_coverage:
        print('{}\t{}'.format(coverage[0], coverage[1]))

    majorlocator = MultipleLocator(1000)
    minorlocator = MultipleLocator(100)
    majorformatter = FormatStrFormatter('%d')

    fig, ax = plt.subplots(1, 1, figsize=(15, 3))
    pos = [i + 1 for i in range(0, len(seq))]
    ticks = [i for i in range(3000, 8000) if i % 100 == 0]

    ax.fill(pos, seq, linewidth=0.75)
    # plt.yscale('log')
    plt.xticks(ticks)
    plt.xlim(3000, 8000)
    plt.ylim(0, 150)

    ax.xaxis.set_major_locator(majorlocator)
    ax.xaxis.set_major_formatter(majorformatter)
    ax.xaxis.set_minor_locator(minorlocator)
    ax.set(xlabel='position', ylabel='Read Count', title='Read Distribution')
    ax.grid()

    # add exon locations
    exon = [('Chr4', 'CDS', 4127, 4149),
            ('Chr4', 'CDS', 4227, 4438),
            ('Chr4', 'CDS', 4545, 4749),
            ('Chr4', 'CDS', 4839, 4901),
            ('Chr4', 'CDS', 4977, 5119),
            ('Chr4', 'CDS', 5406, 5588),
            ('Chr4', 'CDS', 5657, 5855),
            ('Chr4', 'CDS', 6605, 6676),
            ('Chr4', 'CDS', 6760, 6871),
            ('Chr4', 'CDS', 6975, 7056),
            ('Chr4', 'CDS', 7144, 7194),
            ('Chr4', 'CDS', 7294, 7375),
            ('Chr4', 'CDS', 7453, 7638),
            ('Chr4', 'CDS', 7712, 7813),
            ('Chr4', 'CDS', 7914, 7947)]

    span = 8000 - 3000
    for e in exon:
        begin = (e[2] - 3000)/span
        end = (e[3] - 3000)/span
        plt.axhline(5.0, begin, end, color='black', linewidth=6.0)

    plt.show()

    exit(0)

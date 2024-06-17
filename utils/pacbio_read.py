"""=================================================================================================
analyze lengths of sequences in PacBio data, optionally appy a minimum lenght threshold

Michael Gribskov     17 June 2024
================================================================================================="""
import sys
import matplotlib.pylab as plt


def get_data(fasta):
    """---------------------------------------------------------------------------------------------
    read the sequences from the fasta file and stored as
    [ {id:string, length:int} ... ]
    do not store sequences because it may be too big to store in memory

    :param fasta: file      open for reading
    :return: list of dict, int, int     sequence len, minimum and maximum length
    ---------------------------------------------------------------------------------------------"""
    data = []
    len_max = 0
    len_min = 1e100
    sequence = ''
    for line in fasta:
        if line.startswith('>'):
            # title line
            if sequence:
                # only check max and min when th e entry is complete
                len_max = max(len_max, len(sequence))
                len_min = min(len_min, len(sequence))

            entry = {'id': line.rstrip().replace('>', '', 1), 'length': 0}
            data.append(entry)
            sequence = ''
        else:
            # sequence line
            sequence += line.rstrip()
            entry['length'] = len(sequence)

    if sequence:
        # last entry
        len_max = max(len_max, len(sequence))
        len_min = min(len_min, len(sequence))

    return data, len_min, len_max


def get_bins(data, bin_size=100):
    """---------------------------------------------------------------------------------------------
    find the min, max, and average value for each set of bin_size sequences
    
    :param data: list of dict       [{'id', 'length'}, ...]
    :param bin_size: int            number of sequences/bin
    :return: list                   [{'n', 'n_total', 'min', 'ave', 'max', 'sum'}, ...]
    ---------------------------------------------------------------------------------------------"""
    bins = []
    n_total = 0
    for seq in sorted(data, key=lambda l: l['length'], reverse=True):

        if n_total % bin_size:
            # add to current bin, zero means start a new bin
            current['n'] += 1
            current['min'] = min(current['min'], seq['length'])
            current['max'] = max(current['max'], seq['length'])
            current['sum'] += seq['length']
            current['ave'] = current['sum'] / current['n']

        else:
            # start a new bin
            bins.append({
                'n':   1,
                'min': seq['length'],
                'ave': 0,
                'max': seq['length'],
                'sum': seq['length']
                })

            current = bins[-1]

        # end of loop over length data

        n_total += 1
    return bins


def plot(bins):
    """---------------------------------------------------------------------------------------------

    :param bins:
    :return:
    ---------------------------------------------------------------------------------------------"""
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel('Total Length')
    plt.ylabel('Length')
    plt.title(f'Sequence Length Distribution')
    plt.grid(True, linestyle='-', linewidth=0.1)
    h = 0
    ypos = 0
    for b in bins:
        h += 1
        ypos += b['ave'] * b['n']
        plt.vlines(ypos, b['min'], b['max'])

    # plt.text(0.02, .85, 'file: {}\nMean: {:.1f}\nstandard deviation: {:.1f}'.
    #          format(filename, lenmean, lensd), fontsize=10, transform=ax.transAxes)
    # plt.axvline(lenmean, color='red', linewidth=1.5)

    # plt.savefig('{}.pdf'.format(sample))
    plt.show()

    return


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read fasta file, no assumption is made about whether the sequences are on a single line
    fasta = open(sys.argv[1], 'r')
    fasta_data, len_min, len_max = get_data(fasta)
    print(f'sequences read: {len(fasta_data)} ')
    print(f'maximum length: {len_max}')
    print(f'minimum length: {len_min}')

    # collect into bins of fixed number of sequences (instead of fixed bin size)

    bins = get_bins(fasta_data, bin_size=5)
    plot(bins)

    exit(0)

"""=================================================================================================
analyze lengths of sequences in PacBio data, optionally appy a minimum lenght threshold

Michael Gribskov     17 June 2024
================================================================================================="""
import sys


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
                len_max = max( len_max, len(sequence))
                len_min = min( len_min, len(sequence))

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



# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read fasta file, no assumption is made about whether the sequences are on a single line
    fasta = open(sys.argv[1], 'r')
    fasta_data, len_min, len_max = get_data(fasta)
    print(f'{len(fasta_data)} sequences read')

    # collect into bins of fixed number of sequences (instead of fixed bin size)

    bins = get_bins( fasta_data, bin_size = 5 )

    exit(0)

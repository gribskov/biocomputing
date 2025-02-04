"""=================================================================================================
Process a multiple entry fastq file reporting
1) average quality over the entire sequence
2) average quality of just the bases with quality >= q
3) average quality with each sequence truncated at the first base with quality < q
4) the position specific average quality for the untruncated sequence

Not exactly what was implied in the instructions but it succeeds in only traversing the data file
and the quality vectors once

The data structure for the fastq entry is a dictionary with keys {sequence, quality}

Gribskov    4 February 2025
================================================================================================="""


def next_fastq(fp, fastq):
    """---------------------------------------------------------------------------------------------
    return the next sequence in the file in the dit fastq. The quality is the ASCII quality string.
    Newlines are removed from the strings.

    :param fp: file pointer, fastq file open for reading
    :param fastq: dict, two keys: sequence and quality (both strings)
    :return:
    ---------------------------------------------------------------------------------------------"""
    line = fp.readline()  # title, discard for now
    if not line:
        # file is finished
        fastq = {'sequence': None, 'quality': None}
        return False

    sequence = fp.readline().rstrip()
    fp.readline()  # + divide line, discard
    quality = fp.readline().rstrip()

    fastq['sequence'] = sequence
    fastq['quality'] = quality

    return True


def length(fastq):
    """---------------------------------------------------------------------------------------------
    Return the length of the sequence (fastq['sequence'])
    :param fastq: dict      fastq data structure with keys {sequence, quality}
    :return: int            length of sequence
    ---------------------------------------------------------------------------------------------"""
    return len(fastq['sequence'])


def quality(fastq):
    """---------------------------------------------------------------------------------------------
    Convert the ASCII quality string to numeric (Phred quality) and return as a list
    Conversion is q_n = ord(q_ascii) - 33
    
    :param fastq: string    ASCII formatted string, one character per base
    :return: list of int    numerical quality per base position
    ---------------------------------------------------------------------------------------------"""
    qual = []
    for q in fastq['quality']:
        qual.append(ord(q) - 33)

    return qual


def add_counts(quality, stats, qmin):
    """---------------------------------------------------------------------------------------------
    Accumulate counts of number of bases and total quality for all base, high quality bases, and 
    (un) truncated bases. Also accumulate position specific counts and qualities
    
    statistics (stats) are accumulated in a dictionary with keys:
    'nbase_all'    : 0, 'qual_all': 0,       base counts and summed quality, all positions
    'nbase_hq'     : 0, 'qual_hq': 0,        base counts and summed quality, HQ positions only
    'nbase_trunc'  : 0, 'qual_trunc': 0,     base counts and summed quality, truncated at q
    'nbase_by_pos': [], 'qual_by_pos': []    base counts and summed quality, per sequence position

    :param quality: list    numeric quality at each position of sequence
    :param stats: dict      see above
    :param qmin: float      minimum quality to include as HQ, truncation cutoff
    :return: dict           current stats
    ---------------------------------------------------------------------------------------------"""
    # truncated becomes true after the first q < qmin
    truncated = False

    # check to make sure the per position vectors are long enough and extend with zeroes if needed
    if len(stats['nbase_by_pos']) < length(fastq):
        stats['nbase_by_pos'] += [0] * (length(fastq) - len(stats['nbase_by_pos']))
        stats['qual_by_pos'] += [0] * (length(fastq) - len(stats['qual_by_pos']))

    pos = 0
    for q in quality:
        # add values for each position in this fastq entry
        stats['nbase_all'] += 1
        stats['qual_all'] += q
        stats['nbase_by_pos'][pos] += 1
        stats['qual_by_pos'][pos] += q
        pos += 1

        # HQ bases
        if q < qmin:
            # check for truncation
            truncated = True
        else:
            stats['nbase_hq'] += 1
            stats['qual_hq'] += q

        # un-truncated bases
        if not truncated:
            stats['nbase_trunc'] += 1
            stats['qual_trunc'] += q

    return stats


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # minimum quality base for high quality and truncation
    Q_cutoff = 20

    fp = open('8044.5k.fastq', 'r')

    fastq = {}
    stats = {'nbase_all':    0, 'qual_all': 0,
             'nbase_hq':     0, 'qual_hq': 0,
             'nbase_trunc':  0, 'qual_trunc': 0,
             'nbase_by_pos': [], 'qual_by_pos': []
             }

    n_seq = 0
    while next_fastq(fp, fastq):
        # process each fastq entry, accumulating counts in stats
        n_seq += 1
        add_counts(quality(fastq), stats, Q_cutoff)

    print(f'number of sequences: {n_seq}')

    print(f"overall:\t\tbases:{stats['nbase_all']}", end='\t ')
    print(f"average quality: {stats['qual_all'] / stats['nbase_all']:.2f}")

    print(f"high quality:\tbases:{stats['nbase_hq']}", end='\t ')
    print(f"average quality: {stats['qual_hq'] / stats['nbase_hq']:.2f}")

    print(f"truncated:\t\tbases:{stats['nbase_trunc']}", end='\t ')
    print(f"average quality: {stats['qual_trunc'] / stats['nbase_trunc']:.2f}")

    # print in blocks of bsize
    bsize = 20
    print(f'\nPer position quality (blocks of {bsize}):')

    for pos in range(len(stats['qual_by_pos'])):
        end = '\n'
        if (pos + 1) % bsize:
            end = '\t'
        ave = stats['qual_by_pos'][pos] / stats['nbase_by_pos'][pos]
        print(f'{ave:<5.2f}', end=end)

    exit(0)

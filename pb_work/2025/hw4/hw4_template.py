"""
Don't forget the header!
"""


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
    # dummy function  returns 1
    return 1


def quality(fastq):
    # dummy function  returns empty list
    return []


def add_counts(quality, stats, qmin):
    # dummy functions, updates stats dictionary, returns True (what is a more sensible return val)
    return True


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

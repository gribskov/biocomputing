"""=================================================================================================
Homework 3

More practice with dictionaries and functions.

You must use the outline in this template for the general reading of the file.  You may not
change the main program. Your main task is to implement the following functions:

fastq_next:  returns a dictionary containing the next entry in the fastq files.  The keys of the
            dictionary, fq_entry,  are 'id', 'seq', and 'qual'

basecount: calculates the number of each type of base in the fq_entry dictionary and adds i onto
            the dictionary count, optionally applying a minimum threshold for quality

truncate: truncates a fq_entry at the first base with quality lower than threshold and returns a
            dictionary with the same keys but with the sequence and quality truncated

These functions are described in more detail in the docstrings of each function, below

you may create any additional functions you like; the quality function from HW2 can be used
directly.  basecount is essentially the same as composition in HW2.  Truncate is essentially the
same as trim_by_quality.


Michael Gribskov     29 January 2019
================================================================================================="""


def quality(q_string):
    """---------------------------------------------------------------------------------------------
    return a list of quality values based on a quality string.  The quality values are encoded as
    ASCII letters as follows

    quality = ord(letter) - 33

    quality reflects the -10*log10(probability that base call is incorrect) at that position

    :param q_string: string
    :return: list of int
    ---------------------------------------------------------------------------------------------"""
    quality = []
    for q in q_string:
        quality.append(ord(q) - 33)

    return quality


def fastq_next(fq_entry, fq_file):
    """---------------------------------------------------------------------------------------------
    Read the next four lines from the file and returns a dictionary with lines 1, 2, and 4 - the
    ID, sequence, and quality strings, respectively. Trailing newlines are stripped off.

    returns False when an entry cannot be read (i.e., when done), True otherwise

    Dictionary keys: 'id', 'seq', 'qual'

    :param fq_entry: dict
    :param fq_file: filehandle
    :return: bool
    ---------------------------------------------------------------------------------------------"""

    fq_entry['id'] = fq_file.readline().rstrip()
    fq_entry['seq'] = fq_file.readline().rstrip()
    junk = fq_file.readline()
    fq_entry['qual'] = fq_file.readline().rstrip()

    if not fq_entry['qual']:
        return False
    fq_entry['qual'] = quality(fq_entry['qual'])
    return True


def basecount(fq_entry, count, threshold=0):
    """---------------------------------------------------------------------------------------------
    Count the sequence characters, whatever they are, in the sequence.  An additional key 'total'
    holds the total number of bases

    :param fq_entry: dict, keys='id', 'seq', 'qual'
    :param count: dict, keys are DNA bases + 'total'
    :param threshold: int, minimum quality to count
    :return: dict, updated count dict
    ---------------------------------------------------------------------------------------------"""
    if 'total' not in count:
        count['total'] = 0

    seq = fq_entry['seq']
    qual = fq_entry['qual']
    for pos in range(len(seq)):

        if qual[pos] >= threshold:
            count['total'] += 1
            base = seq[pos]
            if base in count:
                count[base] += 1
            else:
                count[base] = 1

    return count


def truncate(fq_entry, threshold):
    """---------------------------------------------------------------------------------------------
    Truncate the sequence and quality at the first position where quality is below threshold

    :param fq_entry: dict, keys='id', 'seq', 'qual'
    :param threshold: int
    :return: dict, keys='id', 'seq', 'qual'
    ---------------------------------------------------------------------------------------------"""
    new_entry = {}

    i = 0
    qual = fq_entry['qual']
    while qual[i] >= threshold:
        i += 1
        if i == len(qual):
            # reached the end of the sequence with good quality
            break

    # i is the index of the first low quality base
    new_entry['id'] = fq_entry['id']
    new_entry['seq'] = fq_entry['seq'][:i]
    new_entry['qual'] = fq_entry['qual'][:i]

    return new_entry


def report(count, title):
    """---------------------------------------------------------------------------------------------
    Print a summary of the base counts and total number of bases

    :param count: dict, keys are bases + 'total'
    :param title: string, printed as a label
    :return: True
    ---------------------------------------------------------------------------------------------"""
    print('\n{}'.format(title))
    for k in count:
        print('\t{}: {}'.format(k, count[k]))

    return True


# --------------------------------------------------------------------------------------------------
# Main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    fastqname = '../HW1/8044.5k.fastq'
    # fastqname = 'test.fq'
    quality_threshold = 20  # pre-adaptation for future command line use

    # Open the file for reading, terminate with unsuccessful status if file cannot be read
    # Defining fastq as None prevents the warning about possibly undefined variable
    fastq = None
    try:
        fq = open(fastqname, 'r')
    except (IOError, OSError):
        print('Error opening file {}'.format(fastqname))
        exit(1)

    all = {}
    high_quality = {}
    trimmed = {}

    fq_entry = {}
    n_entry = 0
    while fastq_next(fq_entry, fq):
        n_entry += 1
        all = basecount(fq_entry, all)
        high_quality = basecount(fq_entry, high_quality, quality_threshold)
        trimmed = basecount(truncate(fq_entry, quality_threshold), trimmed)

    # report
    print('Fastq entries read: {}'.format(n_entry))

    report(all, 'all bases')
    report(high_quality, 'high quality bases, quality >= {}'.format(quality_threshold))
    report(trimmed, 'bases trimmed at quality = {}'.format(quality_threshold))

    exit(0)

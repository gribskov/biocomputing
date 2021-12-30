"""=================================================================================================
mirdeep produces an output file that includes duplicates.  Typically the duplicates are due to
the presence of identical precursor/star/mature miRNA sequences at multiple positions in the
reference genome.  This leads to overcounting of certain RNAs so it is desirable to merge the
identical probes into a single entry

input: all sample CSV file produced by quantify.pl

Michael Gribskov     30 December 2021
================================================================================================="""
import sys


def read_mirdeep_csv(filename):
    """---------------------------------------------------------------------------------------------
    Read the fields into a dict using the column headings in the file as keys.
    Add a field 'type': 'novel' or 'known'

    :param filename:
    :return: list of dict, all columns of information for each miRNA
    ---------------------------------------------------------------------------------------------"""
    fp = None
    try:
        fp = open(filename, 'r')
    except OSError:
        sys.stderr.write(f'read_mirdeep_csv - unable to open file ({filename}')
        exit(1)

    # read down to section titled
    # novel miRNAs predicted by miRDeep2

    for line in fp:
        if line.startswith('novel miRNAs'):
            break

    header = fp.readline()
    title = header.rstrip().split('\t')

    # read the predicted (novel) miRNAs
    mirna = []
    for line in fp:
        if line == '\n':
            # section ends with blank line
            break

        thismirna = {'type':'novel'}
        field = line.rstrip().split('\t')
        f = 0
        for value in field:
            thismirna[title[f]] = field[f]
            f += 1
        mirna.append(thismirna)

    nnovel = len(mirna)
    print(f'Novel miRNAs: {nnovel}')

    # skip to the known miRNA section which begins with
    # mature miRBase miRNAs detected by miRDeep2

    for line in fp:
        if line.startswith('mature miRBase'):
            break

    header = fp.readline()
    title = header.rstrip().split('\t')

    # read the predicted (novel) miRNAs
    for line in fp:
        if line == '\n':
            # section ends with blank line
            break

        thismirna = {'type':'novel'}
        field = line.rstrip().split('\t')
        f = 0
        for value in field:
            thismirna[title[f]] = field[f]
            f += 1
        mirna.append(thismirna)

    nknown = len(mirna) - nnovel
    print(f'Known miRNAs: {nknown}')
    print(f'Total miRNAs: {len(mirna)}')

    return mirna


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    mir = read_mirdeep_csv('../data/result_23_01_2021_t_06_45_02.csv')

    # pre = open('precursor.fa', 'w')
    # mat = open('mature.fa', 'w')
    # star = open('star.fa', 'w')
    #
    # pre.close()
    # mat.close()
    # fp.close()

    exit(0)

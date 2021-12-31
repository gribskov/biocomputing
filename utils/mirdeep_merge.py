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
    Read the fields into a dict using the column headings in the file as keys
    below
    Add a field 'type': 'novel' or 'known'
    Add a field 'sid': with the provisional id or tag id for novel and known miRNAs respectively
    Add a field 'group': used later for skipping merged miRNAs

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

        thismirna = {'type': 'novel', 'group': False}
        field = line.rstrip().split('\t')
        f = 0
        for value in field:
            thismirna[title[f]] = field[f]
            f += 1
        thismirna['sid'] = thismirna['provisional id']
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

    # read the known miRNAs
    for line in fp:
        if line == '\n':
            # section ends with blank line
            break

        thismirna = {'type': 'known', 'group': False}
        field = line.rstrip().split('\t')
        f = 0
        for value in field:
            thismirna[title[f]] = field[f]
            f += 1
        thismirna['sid'] = thismirna['tag id']
        mirna.append(thismirna)

    nknown = len(mirna) - nnovel
    print(f'Known miRNAs: {nknown}')
    print(f'Total miRNAs: {len(mirna)}')

    return mirna


def findmerges(mirna):
    """---------------------------------------------------------------------------------------------
    find which miRNAs are mergable because they have identical precursor sequences

    :param mirna: list of dict, see read_mirdeep_csv()
    :return:
    ---------------------------------------------------------------------------------------------"""
    seqp = ''
    merges = [[-1]]

    for m in sorted(range(len(mirna)), key=lambda x: (mirna[x]['consensus precursor sequence'],
                                                      mirna[x]['sid'])):
        # this sorts by the precursor sequence and the mirdeep id, e.g., 4_6597
        if mirna[m]['consensus precursor sequence'] == seqp:
            merges[-1].append(m)
            mirna[m]['group'] = True

        else:
            if len(merges[-1]) == 1:
                merges[-1] = [m]
            else:
                merges.append([m])
            seqp = mirna[m]['consensus precursor sequence']

    # last item will probably be unmerged
    if len(merges[-1]) == 1:
        merges.pop()

    return merges


def get_id(thismirna):
    """---------------------------------------------------------------------------------------------
    get the primary and alternate ID depending on whether it is a known or novel miRNA

    :param thismirna: dict, one record from mirDeep2 result
    :return: string, string - primary and alternate IDs
    ---------------------------------------------------------------------------------------------"""
    if thismirna['type'] == 'known':
        id = thismirna["mature miRBase miRNA"]
        alt = thismirna['tag id']
    else:
        id = thismirna['provisional id']
        alt = thismirna['example miRBase miRNA with the same seed']

    return id, alt


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    min_mirdeep_score = 0
    min_count = 330

    mirna = read_mirdeep_csv('../data/result_23_01_2021_t_06_45_02.csv')
    merges = findmerges(mirna)

    # print(f'{mirna[-1]}')
    for i in range(len(merges)):
        print(f'\ngroup {i}: {merges[i]}')
        for g in merges[i]:
            m = mirna[g]
            id, alt = get_id(m)

            count = m['total read count']
            seq = m['consensus precursor sequence']
            group = m['group']
            print(f'{g}\t{group}\t{id}\t{alt}\t\t{count}\t\t{seq}')

    pre = open('precursor.fa', 'w')
    mat = open('mature.fa', 'w')
    star = open('star.fa', 'w')

    count_all = {t: 0 for t in ('known', 'novel')}
    count_undup = {t: 0 for t in ('known', 'novel')}
    count_mirdeep = {t: 0 for t in ('known', 'novel')}
    count_count = {t: 0 for t in ('known', 'novel')}
    for m in mirna:

        count_all[m['type']] += 1
        if m['group']:
            continue

        count_undup[m['type']] += 1
        if float(m['miRDeep2 score']) < min_mirdeep_score:
            continue

        count_mirdeep[m['type']] += 1
        if int(m['total read count']) < min_count:
            continue

        count_count[m['type']] += 1

        # write out fasta sequences for those that pass the mirdeep and count cutoffs
        id, alt = get_id(m)
        pre.write(f">{id} {alt}\n{m['consensus precursor sequence']}\n")
        mat.write(f">{id} {alt}\n{m['consensus mature sequence']}\n")
        star.write(f">{id} {alt}\n{m['consensus star sequence']}\n")

    print('!summary')
    print('!Count\tCount')
    print('!Known\tNovel\tID')
    print(f"{count_all['known']}\t\t{count_all['novel']}\t\t total miRNA analyzed")
    print(f"{count_undup['known']}\t\t{count_undup['novel']}\t\t un-duplicated miRNA")
    print(f"{count_mirdeep['known']}\t\t{count_mirdeep['novel']}\t\t mirDeep > {min_mirdeep_score}")
    print(f"{count_count['known']}\t\t{count_count['novel']}\t\t count > {min_count}")

    pre.close()
    mat.close()
    star.close()

    exit(0)

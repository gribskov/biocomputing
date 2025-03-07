"""=================================================================================================
trinity predicted transcripts may be very large due to the presence of many short predicted
transcripts, and transcripts that have few or no counts.

Make a histogram based on salmon counts
cut at a minimum length and merge counts at the desired level: gene, component, or bundle

Michael Gribskov     06 March 2025
================================================================================================="""


def quant_read(infile):
    """---------------------------------------------------------------------------------------------
    read the quant.sf file written by salmon and get the length of each transcript

    :param infile: string       path to quant.sf file (any file will do, all samples are the same)
    :return: dict of dict       dict key is transcript id, inner dict key: length
    ---------------------------------------------------------------------------------------------"""
    quant = open(infile, 'r')
    quant.readline()  # skip column headings

    transcript = {}
    for line in quant:
        if line:
            field = line.split()
            transcript[field[0]] = {'length': float(field[1])}

    quant.close()
    return transcript


def count_read(infile, assembly):
    """---------------------------------------------------------------------------------------------
    read in the merged counts (produced by salmon quantmerge) and add to assembly list

    :param infile: string           path to merged count file
    :param assembly: list of dict   produced by quant_read(), transcript IDs and lengths
    :return:
    ---------------------------------------------------------------------------------------------"""
    count = open(infile, 'r')
    count.readline()  # skip column headings

    count_total = 0
    for line in count:
        field = line.rstrip().split()
        id = field[0]
        tcount = 0
        for c in field[1:]:
            tcount += float(c)

        assembly[id]['count'] = tcount
        count_total += tcount

    count.close()
    return count_total


def merge_count(assembly, level=3):
    """--------------------------------------------------------------------------------------------
    merge counts at the indicated trinity level: 'bundle':1, 'component':2, 'gene': 3, 'isoform':4

    :param assembly: dict of dict   from quant_read() or count_read()
    :param level: int               merging level, default = 3
    :return: dict of dict           key = merged ID, inner dict has list of merged IDs, and count
    --------------------------------------------------------------------------------------------"""
    merged = {}

    for id in assembly:

        # construct the ID at the desired merging level
        field = id.split('_')
        mid = '_'.join(field[:level + 1])
        # print(f'id:{id}\tmerged id:{mid}')

        # add the count into the merged count
        if mid in merged:
            merged[mid]['count'] += assembly[id]['count']
            merged[mid]['member'].append(id)
            merged[mid]['length'] = max(merged[mid]['length'], assembly[id]['length'])
        else:
            merged[mid] = {'member': [id], 'count': assembly[id]['count'], 'length': assembly[id]['length']}

    return merged


def histogram_bylen(total, cut, assembly):
    """---------------------------------------------------------------------------------------------
    calculate the fraction  of reads below each level in cut

    :param total:           total counts, used to set cutoffs for each level
    :param cut: list        cutoff levels for the histogram
    :param assembly: dict   dictionary with transcripts and lengths
    :return: list
    ---------------------------------------------------------------------------------------------"""
    # set up thresholds based on total count
    threshold = []
    lenlevel = []
    countlevel = []
    nlevel = []
    for c in cut:
        threshold.append(c * total)
        lenlevel.append(0)
        countlevel.append(0)
        nlevel.append(0)

    counts = 0
    i = 0
    for t in sorted(assembly, key=lambda l: assembly[l]['length']):
        transcript = assembly[t]
        counts += transcript['count']

        while counts > threshold[i]:
            # find level
            try:
                avct = countlevel[i] / nlevel[i]
            except ZeroDivisionError:
                avct = 0
            print(f'level:{cut[i]}\tcounts:{counts:11.1f} | {countlevel[i]:11.1f}\t'
                  f'n:{nlevel[i]:6d}\tmaxlength:{lenlevel[i]:7.1f}\t avecount:{avct:.1f}')
            i += 1

        nlevel[i] += 1
        countlevel[i] += transcript['count']
        lenlevel[i] = transcript['length']

    try:
        avct = countlevel[i] / nlevel[i]
    except ZeroDivisionError:
        avct = 0

    print(f'level:{cut[i]}\tcounts:{counts:11.1f} | {countlevel[i]:11.1f}\tn:{nlevel[i]:6d}\t'
          f'maxlength:{lenlevel[i]:7.1f}\t avecount:{avct:.1f}\n')

    return lenlevel, countlevel


def histogram_bycount(total, cut, assembly):
    """---------------------------------------------------------------------------------------------
    calculate the fraction  of reads below each level in cut

    :param total:           total counts, used to set cutoffs for each level
    :param cut: list        cutoff levels for the histogram
    :param assembly: dict   dictionary with transcripts and lengths
    :return: list
    ---------------------------------------------------------------------------------------------"""
    # set up thresholds based on total count
    threshold = []
    lenlevel = []
    countlevel = []
    nlevel = []
    ntotal = 0
    for c in cut:
        threshold.append(c * total)
        lenlevel.append(0)
        countlevel.append(0)
        nlevel.append(0)

    counts = 0
    i = 0
    for t in sorted(assembly, key=lambda l: assembly[l]['count']):
        transcript = assembly[t]
        counts += transcript['count']

        while counts > threshold[i]:
            # find level, write out current level when level changes

            try:
                avlen = lenlevel[i] / nlevel[i]
                avct = countlevel[i] / nlevel[i]
            except ZeroDivisionError:
                avlen = 0
                avct = 0

            print(f'level:{cut[i]:<g}\tcounts:{counts:11.1f} | {countlevel[i]:11.1f}\tn:{nlevel[i]:6d}\tavelength:'
                  f'{avlen:7.1f}\t avecount:{avct:.1f}')
            ntotal += nlevel[i]
            i += 1

        nlevel[i] += 1
        countlevel[i] += transcript['count']
        lenlevel[i] += transcript['length']

    try:
        avlen = lenlevel[i] / nlevel[i]
        avct = countlevel[i] / nlevel[i]
    except ZeroDivisionError:
        avlen = 0
        avct = 0

    print(f'level:{cut[i]:<5g}\tcounts:{counts:11.1f} | {countlevel[i]:11.1f}\tn:{nlevel[i]:6d}\tavelength:'
          f'{avlen:7.1f}\t avecount:{avct:.1f}')
    ntotal += nlevel[i]
    # print(f'\ntotal assemblies: {ntotal}\n')

    return lenlevel, countlevel


def filter_bycount(assembly, threshold):
    """---------------------------------------------------------------------------------------------
    keep all assemblies >= threshold

    :param assembly: dict       dictionary with transcripts and lengths
    :param threshold: float     minimum count
    :return: dict               dictionary with selected transcripts and lengths
    ---------------------------------------------------------------------------------------------"""
    selected = {}
    rejected = 0
    for t in sorted(assembly, key=lambda l: assembly[l]['count']):
        transcript = assembly[t]

        if transcript['count'] < threshold:
            rejected += transcript['count']
            continue

        selected[t] = assembly[t]  # not a true copy

    return selected, rejected


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    quantfile = 'data/quant.sf'
    mergefile = 'data/C16C31.numreads.txt'

    # count merge level
    mergelevel = {'bundle': 1, 'component': 2, 'gene': 3, 'isoform': 4}
    merge = mergelevel['gene']
    print(f'Merging at level {merge}')

    # read quant.sf file to get the lengths of the assemblies
    assembly = quant_read(quantfile)
    print(f'assemblies read from {quantfile}: {len(assembly)} ')

    # read the counts and merge with assembly length
    count_total = count_read(mergefile, assembly)
    print(f'total count: {count_total:.1f}')
    merged = merge_count(assembly, merge)
    print(f'merged assemblies: {len(merged)}\n')

    print(f'Length histogram - merged')
    hlen = histogram_bylen(count_total, [0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.02,
                                         0.05, 0.075, 0.10, 0.25, 0.5, 0.75, 0.90, 1.0], merged)

    print(f'Count histogram - merged')
    hcut = histogram_bycount(count_total, [0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.02,
                                           0.05, 0.075, 0.10, 0.25, 0.5, 0.75, 0.90, 1.0], merged)

    # filter by minimum count
    selected, rejected = filter_bycount(merged, 60)
    print(f'\nselected merged transcripts: {len(selected)}  counts removed: {rejected:.1f}'
          f' ({rejected / count_total * 100.0:.1f}%)\n')

    print(f'Length histogram - selected')
    hlen = histogram_bylen(count_total, [0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.02,
                                         0.05, 0.075, 0.10, 0.25, 0.5, 0.75, 0.90, 1.0], selected)

    print(f'Count histogram - selected')
    hcut = histogram_bycount(count_total, [0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.02,
                                           0.05, 0.075, 0.10, 0.25, 0.5, 0.75, 0.90, 1.0], selected)

    # write list of selected transcripts
    selectedfile = 'selected.txt'
    out = open(selectedfile, 'w')
    out.write(f"! transcript                    count        length       members\n")
    for t in selected:
        trans = selected[t]
        out.write(f"!  {t:27s}  {trans['count']:<11.1f}  {trans['length']:<11.1f}  {len(trans['member']):<3d}\n")
        for i in trans['member']:
            out.write(f"{i:30s}  {assembly[i]['count']:<11.1f}  {assembly[i]['length']:<11.1f}\n")

    print(f'Output written to {selectedfile}')
    out.close()

    exit(0)

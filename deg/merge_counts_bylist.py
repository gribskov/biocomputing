"""=================================================================================================


Michael Gribskov     11 March 2025
================================================================================================="""

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    merged_count_file = 'data/C16C31.numreads.txt'
    # selectfile = 'data/count_plus_tax.selected.list'
    selectfile = 'data/count_plus_tax.selected.list'
    outfile = 'data/trinity_selected.merged.count.txt'
    level = 3  # gene level
    merged_count_min = 100

    # read in list of selected transcripts
    select = open(selectfile, 'r')
    selected = []
    sel_n = 0
    for line in select:
        id, rest = line.split(maxsplit=1)
        selected.append(id)
        sel_n += 1

    print(f'selected transcripts from {selectfile}: {sel_n}')
    select.close()

    # read and filter merged read file
    merge = open(merged_count_file, 'r')
    skip_header = merge.readline()
    mlevel = level + 1
    merged = {}
    merged_n = 0
    unmerged_n = 0
    count_total = 0
    rowsum = {}
    for line in merge:
        unmerged_n += 1
        if not unmerged_n%10000:
            print(f'{merged_n} | {unmerged_n}')

        field = line.rstrip().split()
        id = field[0]
        if id not in selected:
            continue

        # merged id at level = level
        token = id.split('_')
        mid = '_'.join(token[:mlevel])
        # print(f'mid:{mid}')

        if mid not in merged:
            merged[mid] = []
            rowsum[mid] = 0
            merged_n += 1
            for col in field[1:]:
                ct = float(col)
                merged[mid].append(ct)
                rowsum[mid] += ct
                count_total += ct

        else:
            # for i in range(1, len(field)):
            for i,ct in enumerate(field[1:]):
                ct = float(ct)
                merged[mid][i-1] += ct
                rowsum[mid] += ct
                count_total += ct

    print(f'\ntotal counts from {merged_count_file}: {count_total:.1f}')
    print(f'transcripts selected: {merged_n}')

    # write out count
    out = open(outfile, 'w')
    out.write(skip_header)
    counts_merged = 0
    merged_n = 0
    for mid in merged:
        if rowsum[mid] >= merged_count_min:
            out.write(f"{mid}")
            for cts in merged[mid]:
                out.write(f'\t{cts:.3f}')
            out.write('\n')
            merged_n += 1

            counts_merged += rowsum[mid]

    print(f'\nmerged counts written to: {outfile}')
    print(f'transcripts above threshold ({merged_count_min}): {merged_n}')
    print(f'merged counts: {counts_merged:.1f}')
    out.close()

exit(0)

"""=================================================================================================
filter a list of trinity IDs, keeping only the ones in a positive list

Michael Gribskov     07 March 2025
================================================================================================="""

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    basefile = 'data/selected.txt'
    posfile = 'data/goodset.out'

    # read positive file into a list of valid IDs, these are the taxonomically valid trinity IDs
    positive = open(posfile, 'r')
    valid_id = []
    for line in positive:
        if line.startswith('!'):
            continue

        field = line.split()
        valid_id.append(field[0])

    print(f'valid IDs: {len(valid_id)} read from {posfile}')
    positive.close()

    # read in list and keep the ID's in the positive list
    selectedfile = 'data/count_plus_tax.selected.list'
    selected = open(selectedfile, 'w')

    rejectedfile = 'data/count_plus_tax.rejected.list'
    rejected = open(rejectedfile, 'w')

    base = open(basefile, 'r')
    select_n = 0
    select_count = 0
    reject_n = 0
    reject_count = 0
    n = 0
    for line in base:
        if line.startswith('!'):
            continue

        n += 1
        # print(f'{n}:{line}', end='')
        if not (n % 1000):
            print('.', end='')
        if not (n % 10000):
            print(f' {n}')

        id, count, length = line.split()
        count = float(count)
        length = float(length)
        if id in valid_id:
            # selected
            selected.write(f'{id}\t{count:.1f}\n')
            select_n += 1
            select_count += count
        else:
            # rejected
            rejected.write(f'{id}\t{count:.1f}\n')
            reject_n += 1
            reject_count += count

    n = select_n + reject_n
    c = select_count + reject_count
    print(f"\nIDs checked: {n} in {basefile}\t{c:11.1f} counts")
    print(
        f'selected: {select_n:11.1f} ({100.0 * select_n / n:5.1f})\tcounts: {select_count:11.1f} ({100.0 * select_count / c:5.1f})')
    print(
        f'rejected: {reject_n:11.1f} ({100.0 * reject_n / n:5.1f})\tcounts: {reject_count:11.1f} ({100.0 * reject_count / c:5.1f})')
    # print(f'rejected: {reject_n:11.1f}\tcounts: {reject_count:11.1f}')

    base.close()
    selected.close()
    rejected.close()

    exit(0)

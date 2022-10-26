"""=================================================================================================
compare to sets of read counts from salmon.  Count files may not have exactly the same genes

Michael Gribskov     26 October 2022
================================================================================================="""
import sys


def open_safe(filename, mode):
    """---------------------------------------------------------------------------------------------
    open a file with check for success. Failure exits with status=1

    :param filename: string     name of file to open
    :param mode: string         mode to open file
    :return: filehandle
    ---------------------------------------------------------------------------------------------"""
    fh = None
    try:
        fh = open(filename, mode)
    except OSError:
        sys.stderr.write(f'compare_counts:opensafe could not open file {filename} in mode {mode}\n')
        exit(1)

    return fh


def read_data(fh):
    """---------------------------------------------------------------------------------------------
    read the multicolumn salmon counts file, return a dictionary with the counts indexed by the
    transcript ID

    :param fh: filehandle   open filehandle
    :return: dict           {transcript id: {column counts}}
    ---------------------------------------------------------------------------------------------"""
    data = {}
    column = fh.readline().rstrip().split()

    for line in fh:
        field = line.rstrip().split()
        transcript = field[0]
        data[transcript] = {}
        for i in range(1, len(field)):
            data[transcript][column[i - 1]] = int(field[i])

    return data


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # open the two files
    f1name = sys.argv[1]
    f2name = sys.argv[2]
    f1 = open_safe(f1name, 'r')
    f2 = open_safe(f2name, 'r')
    print(f'comparing {f1name} and {f2name}')

    f1data = read_data(f1)
    f1len = len(f1data)

    f2data = read_data(f2)
    f2len = len(f2data)

    print(f'counts f1: {f1len}\tf2: {f2len}')

    # find genes in common
    present = {'f1only': [], 'f2only': [], 'both': []}
    allkeys = list(f1data.keys()) + list(f2data.keys())
    allkeys = set(allkeys)
    for t in (allkeys):
        inf1 = t in f1data
        inf2 = t in f2data
        if inf1 and inf2:
            # transcripts in both files
            present['both'].append(t)
        elif inf1:
            present['f1only'].append(t)
        elif inf2:
            present['f2only'].append(t)
        else:
            sys.stderr.write(f'transcript {t} is in neither set')

    print(f"both: {len(present['both'])}\tf1 only: {len(present['f1only'])}\tf2 only: {len(present['f2only'])}")

    f1total_count = {k:0 for k in f1data[present['both'][0]]}
    f2total_count = {k:0 for k in f1data[present['both'][0]]}
    for t in present['both']:
        for c in f1data[t]:
            f1total_count[c] += f1data[t][c]
            f2total_count[c] += f2data[t][c]
            try:
                ratio = f1data[t][c] / f2data[t][c]
            except ZeroDivisionError:
                ratio = 'NA'
            print(f'{t}:{c}\t{f1data[t][c]}\t{f2data[t][c]}\t{ratio}')

    for c in f1total_count:
        print(f'{c}\t{f1total_count[c]}\t{f2total_count[c]}')



    exit(0)

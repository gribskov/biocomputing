"""=================================================================================================
Specifically for annotating trinity counts for viral transcripts in the potato heat shock project
Annotation file has format
!       DN211865_c3_g1  virus   3
TRINITY_DN211865_c3_g1_i1       virus   UniRef90_G8FVS1 1.90e-66        Potato virus S  Coat protein (Fragment) Viruses;Riboviria;Orthorn
TRINITY_DN211865_c3_g1_i1       virus   UniRef90_Q4A3Q8 4.31e-56        Potato virus S  Capsid protein Viruses;Riboviria;Orthorn
TRINITY_DN211865_c3_g1_i1       virus   UniRef90_G8FVR2 8.07e-54        Potato virus S  Coat protein (Fragment) Viruses;Riboviria;Orthorn

count file is from salmon quantmerge
Name    C16T00R1        C16T00R2        C16T00R3        C16T06R1        C16T06R2        C16T06R3        C16T48R1        C16T48R2        C16T48R3        C16T72R1        C16T72R2        C16T72R3        C31T00R1        C31T00R2        C31T00R3        C31T06R1        C31T06R2     C31T06R3 C31T48R1        C31T48R2        C31T48R3        C31T72R1        C31T72R2        C31T72R3
TRINITY_DN89673_c0_g1_i1        13.213  4.459   244.559 154.709 6.019   6.22    91.52   5.982   9.309   81.444  27.307  8.141   223.604 817.639 3393.11 934.259 4434.52 2451.67 9952.48 3221.4  8868.14 4575.28 2140.3  21567.2
TRINITY_DN33670_c2_g1_i1        234.522 98.507  402.017 98.435  49.281  141.128 11.048  334.9   140.119 7       456.503 227.337 96.543  305.948 673.68  147.729 266.053 132.682 256.789 97.972  160.422 149.072 170.393 99.334
TRINITY_DN33670_c0_g1_i1        0       0       0       0       22      197.043 0       0       0       0       359.075 0       309.813 530.507 738.631 327.636 732.952 620.672 1858.51 395.009 1380.24 1063.26 472.044 625.606
TRINITY_DN7062_c0_g2_i1 0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       2       0       0       0       0       0


Michael Gribskov     02 May 2025
================================================================================================="""
import sys
from collections import defaultdict


def annotation_read(fname):
    """---------------------------------------------------------------------------------------------
    read the annotation file
    :param fname: str       path to file with annotation data
    :return: dict           indexed by truncated trinity id
    ---------------------------------------------------------------------------------------------"""
    infile = open(fname, 'r')
    anno = defaultdict(lambda: defaultdict(int))
    for line in infile:
        if line.startswith('!'):
            continue

        field = line.split('\t')
        tid = trinity_truncate(field[0], 3)
        anno[tid][field[4]] += 1

    infile.close()
    return anno


def count_read(fname):
    """---------------------------------------------------------------------------------------------
    read the annotation file

    :param fname: str       path to file with count data
    :return: dict           indexed by truncated trinity id
    ---------------------------------------------------------------------------------------------"""
    infile = open(fname, 'r')

    # get number of columns from header
    field = infile.readline().split()
    ncol = len(field) - 1

    count = defaultdict(lambda: [0 for _ in range(ncol)])
    for line in infile:
        field = line.rstrip().split('\t')
        tid = trinity_truncate(field[0], 3)

        for i in range(ncol):
            count[tid][i] += round(float(field[i + 1]), 1)

    infile.close()
    return count


def trinity_truncate(idstr, level=3):
    """---------------------------------------------------------------------------------------------
    truncate the trinity id at level. Return input string if level == 0
    level   id
    4       TRINITY_DN223755_c0_g1_i1
    3       TRINITY_DN223755_c0_g1
    2       TRINITY_DN223755_c0
    1       TRINITY_DN223755

    :param idstr: str       trinity ID
    :param level: int       truncation level
    :return: str            truncated ID
    ---------------------------------------------------------------------------------------------"""
    if level == 0:
        return idstr

    return '_'.join(idstr.split('_')[:level + 1])


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    anno = annotation_read(sys.argv[1])
    print(f'annotations read from {sys.argv[1]}: {len(anno)}')
    count = count_read(sys.argv[2])
    print(f'counts read from {sys.argv[2]}: {len(anno)}')

    out = open(sys.argv[3], 'w')

    virus = defaultdict(int)
    for tid in anno:
        for v in anno[tid]:
            virus[v] += anno[tid][v]

    print(f'\nCommon virus labels')
    out.write(f'\nCommon virus labels\n')
    n = 0
    ascii = 97
    common = defaultdict(lambda: '')
    for v in sorted(virus, key=lambda v: virus[v], reverse=True):
    # print(f'{v}: {virus[v]}')
        common[v] = chr(ascii)
        n += 1
        if n > 6:
            ascii = 111
        else:
            print(f'{common[v]}\t {v}')
            out.write(f'{common[v]}\t {v}\n')
        ascii += 1

    print()
    out.write('\n')

    rowsum = defaultdict(int)
    for tid in count:
        rowsum[tid] = sum(count[tid])

    for tid in sorted(rowsum, key=lambda k: rowsum[k], reverse=True):
    # construct an abbreviated annotation string using the common virus
    # abbreviations from above
        vstr = ''
        for v in sorted(anno[tid], key=lambda x: anno[tid][x], reverse=True):
            vstr += f'{common[v]}{anno[tid][v]}'

        values = ''
        for v in count[tid]:
            values += f'\t{round(v, 1)}'

        print(f'{tid}\t{vstr}\t{round(rowsum[tid], 1)}\t{values}')
        out.write(f'{tid}\t{vstr}\t{round(rowsum[tid], 1)}\t{values}\n')

    out.close()
    exit(0)

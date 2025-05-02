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


def annotation_read():
    """---------------------------------------------------------------------------------------------
    read the annotation file
    :return:
    ---------------------------------------------------------------------------------------------"""
    infile = open(sys.argv[1], 'r')
    anno = defaultdict(lambda: defaultdict(int))
    for line in infile:
        if line.startswith('!'):
            continue

        field = line.split('\t')
        tid = trinity_truncate(field[0], 3)
        anno[tid][field[4]] += 1

    infile.close()
    return anno


def count_read():
    """---------------------------------------------------------------------------------------------
    read the annotation file
    :return:
    ---------------------------------------------------------------------------------------------"""
    pass


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
    anno = annotation_read()
    count = count_read()

    exit(0)

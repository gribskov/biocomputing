"""=================================================================================================
read the npinter experimentally validated datafile
http://bigdata.ibp.ac.cn/npinter5

Michael Gribskov     08 September 2025
================================================================================================="""
import os
# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # path = 'A:/mrg/npinter/interaction_NPInterv5.expr.txt'
    path = 'A:/mrg/npinter/lncRNA_interaction.txt'

    rna = open(path, 'r', encoding='utf-8')

    itype = {}
    mlist = {}
    n = 0
    nn = 0
    for line in rna:
        n += 1
        field = line.split('\t')


        if field[14] != 'RNA-Protein':
            continue

        print(line)
        i = 0
        for f in field:
            print(f'\t{i}:{f}')
            i += 1

        nn += 1
        inter = field[14]
        method = field[15]
        print(method)
        if inter in itype:
            itype[inter] += 1
        else:
            itype[inter] = 1
        if method in mlist:
            mlist[method] += 1
        else:
            mlist[method] = 1

    print(n,nn)
    for i in itype:
        print(f'{i}: {itype[i]}')

    print(f'\nmethods')
    for i in mlist:
        print(f'{i}: {mlist[i]}')

    exit(0)

"""=================================================================================================
gathers the results of the blast search to list the annotations by trinity gene or component
================================================================================================="""
import sys

info = {}
ngene = 0
for line in sys.stdin:

    # remove everything after n=
    nspot = line.find(' n=')
    # split the first field
    field = line[:nspot].lower().split()
    trinity = field[0].split('_')
    # gene level gene = '{}_{}_{}'.format(trinity[1], trinity[2], trinity[3] )
    # component
    gene = '{}_{}'.format(trinity[1], trinity[2] )

    lastfield = len(field) - 1
    if field[lastfield].find('fragment') > -1:
        lastfield -= 1
    if field[lastfield].find('protein') > -1:
        lastfield -= 1

    descrip = ' '.join(field[12:lastfield+1])
    #print(gene, descrip)
    if gene in info:
        if descrip not in info[gene]:
            info[gene].append(descrip)
    else:
        info[gene] =  [descrip]
        ngene += 1

print('{} genes'.format(ngene))
for gene in info:
    print('{}\t{}'.format(gene, info[gene][0]))
    for i in range(1,len(info[gene])):
        print('\t{}'.format(info[gene][i]))

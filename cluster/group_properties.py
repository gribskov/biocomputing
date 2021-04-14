"""=================================================================================================
Using the information in the fasta file, annotate the groups produced by blast_singlelinkage.py

Michael Gribskov     14 April 2021
================================================================================================="""
import sys


def readcluster(cluster):
    """---------------------------------------------------------------------------------------------
    Read cluster file.  Sample input
    group 8997	19 members:
	    Blessica_CDS_106, Catdawg_CDS_109, Corndog_CDS_105, Familton_CDS_109, Firecracker_CDS_108...
    group 8998	2 members:
        Blessica_CDS_107, Familton_CDS_110
    ...
    36470727 hits clustered
    314592 sequences
    groups size histogram (number of sequences, number of groups
    1	15580
    2	4664
    3	2256
    ...

    :param cluster: filehandle
    :return:
    ---------------------------------------------------------------------------------------------"""
    group_list = []
    for line in cluster:
        if line.endswith('hits clustered\n'):
            break

        if line.startswith('group'):
            field = line.rstrip().split()
            group = int(field[1])
            group_n = int(field[2])
        else:
            names = line.strip().split(', ')
            namecount = {}
            for gene in names:
                phage, genename = gene.split('_', maxsplit=1)
                if genename in namecount:
                    namecount[genename] += 1
                else:
                    namecount[genename] = 1

            group_list.append({'name':group, 'n':group_n, 'member':names, 'gene':namecount})

    return group_list


def read_doc(title):
    """-----------------------------------------------------------------------------------------
    from the fasta file titles, store the name of the gene and the documentation as a dict

    >20ES_CDS_10 lysin B
    >20ES_CDS_11 terminase
    >20ES_CDS_12 portal protein
    >20ES_CDS_13 capsid maturation protease
    >20ES_CDS_14 scaffold protein

    :param title: filehandle, list of title of sequences from grep '>' file.fa
    :return: dict, phage_name:doc
    -----------------------------------------------------------------------------------------"""
    doc = {}
    for line in title:

        name, function = line.strip('>\n').split(' ', maxsplit=1)
        if not function:
            function = 'unk'

        doc[name] = function

    return doc


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    title = open(sys.argv[1], 'r')
    cluster = open(sys.argv[2], 'r')
    out = open(sys.argv[3], 'w')

    groups = readcluster(cluster)
    doc = read_doc(title)

    function = {}
    gnum = 0
    for group in groups:
        fxn = {}
        for member in group['member']:
            f = doc[member]
            if f in fxn:
                fxn[f] += 1
            else:
                fxn[f] = 1

            if f in function:
                function[f]['count'] += 1
            else:
                function[f] = {'count':1, 'groupnum':[]}

        group['function'] = {f:fxn[f] for f in sorted(fxn, key=lambda x: fxn[x], reverse=True)}


        # store the group number in the list of groups with this function
        # because some groups are merged, the group number in the list does not match group['name']
        for f in fxn:
            function[f]['groupnum'].append(gnum)

        gnum += 1

    # summarize functions
    out.write('groups by function\n\n')
    for fxn in sorted(function):
        out.write('function:{}\tsequences:{}\n'.format(fxn, function[fxn]['count']))
        for g in function[fxn]['groupnum']:
            otherfxn = str(groups[g]['function']).replace("'", "").replace(': ', ':').strip('{}')
            out.write('\tgroup:{:<6d}\tsequences:{:<6d}\tthis:{:<6d}\tall:{}\n'.format(
                groups[g]['name'],groups[g]['n'], groups[g]['function'][fxn], otherfxn))

    # summarize groups
    out.write('\nfunctions by groups\n\n')
    for group in sorted(groups, key=lambda g: g['n'], reverse=True):
        out.write('group:{}\tsequences:{}\n'.format(group['name'], group['n']))
        fxn_num = []
        for f in sorted(group['function'], key=lambda f: group['function'][f], reverse=True):
            fxn_num.append('{}={}'.format(f, group['function'][f]))

        out.write('\tfunction:{}\n'.format('; '.join(fxn_num)))

        gene_num = []
        for g in sorted(group['gene'], key=lambda f: group['gene'][f], reverse=True):
            gene_num.append('{}={}'.format(g, group['gene'][g]))

        out.write('\tgene:{}\n'.format('; '.join(gene_num)))


    exit(0)

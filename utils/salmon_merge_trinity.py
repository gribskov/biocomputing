"""=================================================================================================
merge salmon output at the cluster, component, and gene, and isoform levels
works on multiple files to produce a counts file for DESeq2
================================================================================================="""
import sys
import os
import glob
import re

# regex for splitting Trinity IDs
idre = re.compile(r'>*TRINITY_([^_]+)_c(\d+)_g(\d+)_i(\d+)')


def id_split(id):
    """---------------------------------------------------------------------------------------------
    split the Trinity ID into its parts and return
    the string "TRINITY_" is removed

    :param id: str
    :return: cluster, component, group, isoform
    ---------------------------------------------------------------------------------------------"""
    match = idre.match(id)

    return match.groups()


def combine(id, data, count):
    """---------------------------------------------------------------------------------------------
    Add the count to the data dict
    :param id: str
    :param data: dict, list of counts
    :param count: salmon count
    :return: number of id's in list
    ---------------------------------------------------------------------------------------------"""
    if id in data:
        data[id] += count
    else:
        data[id] = count

    return len(data)


# main =============================================================================================

if __name__ == '__main__':

    pathtarget = './'
    filetarget = 'salmon.all.*'
    target = pathtarget + '/' + filetarget
    print('target:', target)

    for filename in glob.glob(target):

        try:
            salmon = open(filename, 'r')
        except IOError:
            sys.stderr.write('unable to open salmon file ({})\n'.format(filename))
            exit(1)

        # get the column name from the directory name.  the assumption is the directory is named
        # something on_merge_trin   path = os.path.split(filename)
        # dir = os.path.split(path[0])
        # column = dir[-1]
        column = os.path.split(os.path.split(filename)[0])[-1]
        column = column.replace('salmon.', '')
        sys.stdout.write('\nfile: {}\tcolumn: {}\n'.format(filename, column))

        salmon = open(filename, 'r')
        # title line
        salmon.readline()

        # read and store the data
        isoform = {}
        gene = {}
        component = {}
        cluster = {}
        n_isoform = n_gene = n_component = n_cluster = 0
        for line in salmon:
            n_isoform += 1
            (raw_id, effective, length, TPM, count) = line.rstrip().split()
            count = float(count)

            # generate names at each hierarchical level
            (cl, co, g, i) = id_split(raw_id)
            isoform_id = '{}_c{}_g{}_i{}'.format(cl, co, g, i)
            gene_id = '{}_c{}_g{}'.format(cl, co, g)
            component_id = '{}_c{}'.format(cl, co)
            cluster_id = '{}'.format(cl)

            isoform[isoform_id] = count
            n_gene = combine(gene_id, gene, count)
            n_component = combine(component_id, component, count)
            n_cluster = combine(cluster_id, cluster, count)

        sys.stdout.write('\t{} isoforms\n'.format(n_isoform))
        sys.stdout.write('\t{} genes\n'.format(n_gene))
        sys.stdout.write('\t{} components\n'.format(n_component))
        sys.stdout.write('\t{} clusters\n'.format(n_cluster))

    exit(0)

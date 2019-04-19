"""=================================================================================================
merge salmon output at the cluster, component, and gene, and isoform levels
works on multiple files to produce a counts file for DESeq2
================================================================================================="""
import sys
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

    filename = '../data/salmon.quant.sf'
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

    sys.stdout.write('{} isoforms read from {}\n'.format(n_isoform, filename))
    sys.stdout.write('{} genes\n'.format(n_gene))
    sys.stdout.write('{} components\n'.format(n_component))
    sys.stdout.write('{} clusters\n'.format(n_cluster))

    exit(0)

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


def read_and_split(salmon_file_name):
    """---------------------------------------------------------------------------------------------
    Read the salmon output and return dicts at each hierarchical level.

    :param salmon: filehandle
    :return: dicts: cluster, component, gene, isoform
    ---------------------------------------------------------------------------------------------"""
    try:
        salmon = open(salmon_file_name, 'r')
    except IOError:
        sys.stderr.write('read_and_split - unable to open salmon file ({})\n'.format(
            salmon_file_name))
        exit(1)

    # skip the column header
    salmon.readline()

    # read and store the data
    isoform = {}
    gene = {}
    component = {}
    cluster = {}
    n_isoform = n_gene = n_component = n_cluster = 0
    for line in salmon:
        n_isoform += 1
        # if n_isoform > 1000:
        #     break

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

    salmon.close()
    return cluster, component, gene, isoform


def write_counts(all, prefix, threshold):
    """---------------------------------------------------------------------------------------------
    Apply threshold and writeout merged counts

    :param all: dict of dicts, each dict is one hierarchical level
    :param threshold:
    :return: dict of counts at each level
    ---------------------------------------------------------------------------------------------"""

    prefix = "merge_salmon_"
    outstr = ''
    threshold = 50

    # couns of transcripts that pass threshold
    passed = {'isoform': 0, 'gene': 0, 'component': 0, 'cluster': 0}

    for level in ['isoform', 'gene', 'component', 'cluster']:
        fname = '{}{}.tsv'.format(prefix, level)
        output = open(fname, 'w')

        # column headings
        for column in sorted(all):
            output.write('\t{}'.format(column))
        output.write('\n')

        for id in all[column_name][level]:
            outstr = '{}\t'.format(id)
            sum = 0
            for col in sorted(all):
                outstr += '\t{}'.format(round(all[col][level][id]))
                sum += all[col][level][id]

            if sum > threshold:
                passed[level] += 1
                output.write('{}\n'.format(outstr))
        output.close()

    return passed


# main =============================================================================================

if __name__ == '__main__':

    pathtarget = '/scratch/snyder/m/mgribsko/avocado/trinity/salmon/*.salmon'
    filetarget = 'quant.sf'
    target = pathtarget + '/' + filetarget

    all = {}
    for filename in glob.glob(target):

        # get the column name from the directory name.  the assumption is the directory is named
        # something like sample.salmon
        column_name = os.path.split(os.path.split(filename)[0])[-1]
        sys.stdout.write('\nfile: {}\tcolumn: {}\n'.format(filename, column_name))

        # read in all counts and store in all{}

        cluster, component, gene, isoform = read_and_split(filename)
        all[column_name] = {'isoform': isoform,
                            'gene': gene,
                            'component': component,
                            'cluster': cluster}

        # report for this file
        sys.stdout.write('\t{} isoforms\n'.format(len(isoform)))
        sys.stdout.write('\t{} genes\n'.format(len(gene)))
        sys.stdout.write('\t{} components\n'.format(len(component)))
        sys.stdout.write('\t{} clusters\n'.format(len(cluster)))

        # threshold counts and write to files
        # the sum of counts must be > threshold
        prefix = "merge_salmon_"
        threshold = 50
        passed = write_counts(all, prefix, threshold)

        all[column_name] = {'isoform': isoform, 'gene': gene, 'component': component,
                            'cluster': cluster}
        # end of loop over count files

    sys.stdout.write('\npassed threhold={}\n'.format(threshold))
    for level in ['isoform', 'gene', 'component', 'cluster']:
        sys.stdout.write('\t{1}\t{0}s\n'.format(level, passed[level]))


    exit(0)

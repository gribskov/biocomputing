"""=================================================================================================
Construct the transcript-to-gene map for use with --transcript-to-gene-map option
input: trinity predicted transcripts
level: level to cluster transcripts
    cluster (e.g. TTRINITY_DN88461)
    component (e.g. RINITY_DN88461_c0)
    gene (e.g. TRINITY_DN88461_c0_g1)
    isoform (transcript to gene map not needed)

for the transcript ID: TRINITY_DN88461_c0_g1_i1

Usage:
    rsem_gene_map_trinity.py <transcript_file> <level>


Michael Gribskov     01 June 2018
================================================================================================="""
import sys


def add(target, key, value):
    """---------------------------------------------------------------------------------------------
    add 1 to the count of a key in a dict, target.  Create key if unknown
    :param target: dictionary
    :param key: string, dictionary key
    :param value: string, complete id
    :return: True
    ---------------------------------------------------------------------------------------------"""
    if key in target:
        target[key].append(value)
    else:
        target[key] = [value]

    return True


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    trinity_file = sys.argv[1]
    level = sys.argv[2]

    sys.stderr.write('rsem_gene_map_trinity\n\tinput: {}\n\tlevel: {}\n\n'.
                     format(trinity_file, level))

    try:
        trinity = open(trinity_file, 'r')
    except:
        sys.stderr.write('Error opening trinity file ({})'.format(trinity_file))

    clusterlist = {}
    componentlist = {}
    genelist = {}
    isoformlist = {}
    nline = 0
    nisoform = 0

    field = []

    for line in trinity:
        nline += 1
        if not line.startswith('>'):
            continue

        nisoform += 1
        id, info = line.replace('>', '').split(maxsplit=1)
        field = id.split('_')
        cluster = '{}_{}'.format(field[0], field[1])
        component = '{}_{}'.format(cluster, field[2])
        gene = '{}_{}'.format(component, field[3])
        isoform = '{}_{}'.format(gene, field[4])

        add(clusterlist, cluster, isoform)
        add(componentlist, component, isoform)
        add(genelist, gene, isoform)

    sys.stderr.write('lines read: {}\n'.format(nline))
    sys.stderr.write('clusters: {}\n'.format(len(clusterlist)))
    sys.stderr.write('components: {}\n'.format(len(componentlist)))
    sys.stderr.write('genes: {}\n'.format(len(genelist)))
    sys.stderr.write('isoforms: {}\n'.format(nisoform))

    if level == 'cluster':
        group = clusterlist
    elif level == 'component':
        group = componentlist
    elif level == 'gene':
        group = genelist
    else:
        sys.sterr.write('unknown level: {}'.format(level))

    for id in group:
        for mapped in group[id]:
            sys.stdout.write('{}\t{}\n'.format(id, mapped))

    exit(0)

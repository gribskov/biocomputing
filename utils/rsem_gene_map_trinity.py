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

    field = []

    for line in trinity:
        nline += 1
        if not line.startswith('>'):
            continue

        id, info = line.replace('>', '').split(maxsplit=1)
        field = id.split('_')
        cluster = '{}_{}'.format(field[0], field[1])
        component = '{}_{}'.format(cluster, field[2])
        gene = '{}_{}'.format(component, field[3])

        print('{}\n\t{}\n\t{}\n\t{}\n'.format(id, cluster, component, gene))

    sys.stderr.write('lines read: {}'.format(nline))

    exit(0)

"""=================================================================================================
construct trinity gene map at the gene, component, or cluster level

TRINITY_DN82208_c0_g1_i1
--------cluster         all components, genes, and isoforms from the same trinity cluster (TRINITY_DN82208)
---------component      all genes, and isoforms with same cluster and component (TRINITY_DN82208_c0)
----------------gene    all isoforms with same cluster, component, and gene (TRINITY_DN82208_c0_g1)

if you want isoforms you don't really need the gene  map

usage:
    trinity_genemap.py  Trinity.fasta t|c|g
================================================================================================="""
import sys

# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    usage = 'trinity_genemap.py  Trinity.fasta t|c|g'
    trinity_file = sys.argv[1]
    level = sys.argv[2]

    try:
        trinity = open(trinity_file, 'r')
    except:
        sys.stderr.write('Error opening input file ({}\n)'.format(trinity_file))
        sys.stderr.write(usage + '\n')

    # TRINITY_DN82208_c0_g1_i1
    col = 2
    if level == 'c':
        col = 3
    elif level == 'g':
        col = 4
    elif level == 'i':
        col = 5

    for line in trinity:
        if line.startswith('>'):
            line = line[1:]
            title, dummy = line.split(' ', maxsplit=1)
            field = title.split('_')

            sys.stdout.write('{}\t{}\n'.format('_'.join(field[:col]),'_'.join(field[:5])))


exit(0)

"""==============================================================================================
given a databases search (blastx) with a trinity assembly with IDs like DN88418_c8_g2_i2
for each cluster (first level trinity group, e.g. DN88418
get the list of hits that map to any member of the cluster
for each pair of transcripts measure the distance over each possible match
    residues that overlap in both
    or residues that match in the same protein (differences may e indicate isoforms or fragments)
cluster to identify disjoint sets
    hierarchical, measure distance between and within groups
    kmeans ? how to do it with distances
        test for best number of groups
================================================================================================="""
import sys

def distance(cluster):
    """---------------------------------------------------------------------------------------------
    calculate distance between members of the cluster

    :param cluster:
    :return:
    ---------------------------------------------------------------------------------------------"""
    return True

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    blastfile = sys.argv[1]
    try:
        blast = open(blastfile, 'r')
    except OSError:
        sys.stderr.write('Unable to open blast search file ({}'.format(blastfile))
        exit(1)

    ecutoff = 1e-40
    clusters = {}
    ncluster = 0
    nline = 0
    for line in blast:
        line = line.rstrip()
        nline += 1
        piece = line.split()
        tag, cluster, component, gene, isoform = piece[0].split('_')
        info = {'id': piece[0],
                'qlen':int(piece[1]),
                'qbegin': int(piece[2]),
                'qend': int(piece[3]),
                'sid': piece[4],
                'slen': int(piece[5]),
                'sbegin': int(piece[6]),
                'send': int(piece[7]),
                'evalue':float(piece[10]),
                'descrip': ' '.join(piece[12:])
                }
        if info['evalue'] > ecutoff:
            continue

        if cluster in clusters:
            clusters[cluster].append(info)
        else:
            clusters[cluster] = [info]
            ncluster += 1

        if ncluster > 1:
            break

    for cluster in clusters:
        dist = distance(cluster)

    sys.stdout.write('{} lines read from {}'.format(nline, blastfile))

    exit(0)
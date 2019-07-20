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

blastfile = sys.argv[1]
try:
    blast = open(blastfile, 'r')
except OSError:
    sys.stderr.write('Unable to open blast search file ({}'.format(blastfile))
    exit(1)

nline = 0
for line in blast:
    nline += 1
    piece = line.split()

sys.stdout.write('{} lines read from {}'.format(nline, blastfile))
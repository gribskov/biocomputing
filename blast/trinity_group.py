'''
trinity_group.py
based on the results of a blast search choose the best isoform of each assembly

based on Trinity names, e.g. 
TRINITY_DN58810_c2_g3_i3
'''
import argparse
import re

# defaults

commandline = argparse.ArgumentParser( 
    description='choose the best Trinity assemblies based on results of diamond blastx search ')
commandline.add_argument( 'trinity_idx',
                          help='index for trinity names',
                          type=argparse.FileType('r'))
commandline.add_argument( 'blast',
                          help='tabular blast search result',
                          type=argparse.FileType('r'))

# process command line arguments
cl = commandline.parse_args()

print( '\ntrinity_group.py - choose best isoforms\n' )

n_trinity = 0
trinity_list = {}
for line in cl.trinity_idx:
    n_trinity += 1
    left, right = line.rstrip().split()
    trinity_list[ right ] = left.lstrip('>')

print( n_trinity, 'isoforms read from', cl.trinity_idx.name )
cl.trinity_idx.close()

#for tname in trinity_list:
#    print( '     left:', tname )
#    print( '     right:', trinity_list[tname] )

# s572    UPI000981EDA0   83.5    67.9    1750    6.1e-192
# qseqid  sseqid          pident  qcovhsp score   evalue
n_blast = 0;
matches = {}
for line in cl.blast:
    line = line.rstrip()
    #print( line )
    n_blast += 1
    qseqid, sseqid, pident, qcov, score, evalue = line.split()
    isoform = trinity_list[qseqid]
    if not isoform in matches:
        matches[isoform] = []
    matches[isoform].append({'sseqid':sseqid, 'pident':pident, 'qcov':qcov, 'score':score, 'evalue':evalue})

print( n_blast, 'blast results read from', cl.blast.name )

for qid in matches:
    print(qid)
    for hit in matches[qid]:
        for key in hit:
            print('    ',key,':', hit[key],end='')

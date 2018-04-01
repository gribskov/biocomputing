"""=================================================================================================
Plot the distribution of reads on a reference sequence based on mapped reads in a SAM file

SAM format is (all one line, whitespace separated fields)

SRR5295840.120  163     AT1G07250.1     101     44      1S150M  =       347     398     NCT...CAG
#A<...FJF AS:i:300     XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:150        YS:i:300        YT:Z:CP

SRR5295840.120  QNAME   read name
163             FLAG    mapping bit flags
AT1G07250.1     RNAME   reference sequence name
101             POS     leftmost position of mapped read
44              MAPQ    mapping quality
1S150M          CIGAR   alignment
=               RNEXT   name of mate/next read
347             PNEXT   position of mate/next read
398             TLEN    inferred insert size
NCT...CAG       SEQ     sequence
#A<...FJF       QUAL    quality

the remaining columns are application specific
AS:i:300     XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:150        YS:i:300        YT:Z:CP
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

map = None
try:
    map = open(sys.argv[1], 'r')
except:
    print('unable to open input file ({}'.format(sys.argv[1]))
    exit(1)

nread = 0
seq = []

# assume that the mapped read begins at POS
# use the CIGAR string to increment counts in the seq array

for line in map:
    nread += 1
    field = line.split()
    # print('{}\t{}'.format(field[0], field[8]))
    pos = int(field[3])
    cigar = field[5]

    if nread > 1000:
        break
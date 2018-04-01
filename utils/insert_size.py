"""=================================================================================================
insert_size.py

Calulate insert size based on a SAM file of mapped reads. To get only high quality mapped read pairs
use the samtools command samtools view -q 20 -f 0x82 SRR5295840.bam > SRR5295840.mapped
SAM format is (all one line, whitespace separated fields)

SRR5295840.120  163     AT1G07250.1     101     44      1S150M  =       347     398     NCT...CAG
#A<...FJF AS:i:300     XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:150        YS:i:300        YT:Z:CP

the insert length is field 8, the last field before the sequence

Michael Gribskov     1 April 2018
================================================================================================="""
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import statistics as stat

map = None
try:
    map = open(sys.argv[1], 'r')
except:
    print('unable to open input file ({}'.format(sys.argv[1]))
    exit(1)

nread = 0
lendata = []
for line in map:
    nread += 1
    field = line.split()
    # print('{}\t{}'.format(field[0], field[8]))
    insert = float(field[8])
    insert = max( insert, -insert)
    lendata.append(insert)

    if nread > 100000:
        break

lenmean = stat.mean(lendata)
lensd = stat.stdev(lendata)

print('\n{} reads read from {}'.format(nread, sys.argv[1]))
print('mean: {:.3f}\tstandard devation: {:.3f}'.format( lenmean, lensd))

fig = plt.figure()
ax = fig.add_subplot(111)

n, bins, patches = plt.hist(lendata, bins=100, normed=1,
                            facecolor='blue', edgecolor='black', linewidth=0.25, alpha=0.75 )

plt.xlabel('Length')
plt.ylabel('Probability')
plt.title('Library Insert Length')
plt.grid(True, linestyle='-', linewidth=0.1)

plt.text( 0.02, 0.9, 'mean: {:.1f}\nstandard deviation: {:.1f}'.format(lenmean, lensd), fontsize=10, transform=ax.transAxes)
plt.axvline( lenmean, color='red', linewidth=1.5)

plt.show()

exit(0)
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
import os.path as path
import re
import matplotlib.pylab as plt
import numpy as np
from scipy import stats as stat

# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':
    # TRIMFRAC = 0.0025
    TRIMFRAC = 0.005
    MAXREADS = 1000000 - 1
    MAXREADLEN = 1200
    BINSIZE = 20

    map = None
    try:
        map = open(sys.argv[1], 'r')
    except:
        print('unable to open input file ({}'.format(sys.argv[1]))
        exit(1)

    # trim filename
    filename = path.basename(sys.argv[1])
    sample = re.sub('\.[\w]+$', '', filename)
    print('Calculating average inset size for: {}'.format(sample))

    nread = 0
    lendata = []
    interval = int(MAXREADS / 20)
    hist = [0 for _ in range(MAXREADLEN+1)]
    for line in map:
        if line.startswith('@'):
            # skip header lines
            continue

        field = line.split()
        insert = float(field[8])
        insert = max(insert, -insert)       # invert negative values
        if insert == 0 or insert > MAXREADLEN:
            continue

        nread += 1
        if not nread % interval:
            # screen trace
            print('.', end=' ', flush=True)


        lendata.append(insert)
        hist[int(insert)] += 1

        if nread > MAXREADS:
            break
    print()

    tdata = stat.trimboth(lendata, TRIMFRAC)
    descrip = stat.describe(tdata)
    lenmean = descrip.mean
    lensd = np.sqrt(descrip.variance)

    print('\n{} reads read from {}'.format(nread, filename))
    print('\n{} fraction of values trimmed from each end'.format(TRIMFRAC))
    print('trimmed mean: {:.3f}\ttrimmed standard devation: {:.3f}'.format( lenmean, lensd))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # insert size histogram
    bins = [i for i in range(int(descrip.minmax[1])) if not i % BINSIZE]
    bins.append(int(descrip.minmax[1]))
    n, bins, patches = plt.hist(tdata, normed=True, bins=bins,
                                facecolor='blue', edgecolor='black', linewidth=0.25, alpha=0.75 )
    # Gaussian fit
    xt = plt.xticks()[0]
    xmin = min(xt)
    xmax = max(xt)
    lnspc = np.linspace(xmin, xmax, len(tdata))
    pdf = stat.norm.pdf(lnspc, lenmean, lensd)
    plt.plot(lnspc, pdf)

    plt.xlabel('Length')
    plt.ylabel('Probability')
    plt.title('{} Insert Length'.format(sample))
    plt.grid(True, linestyle='-', linewidth=0.1)

    plt.text( 0.02, .85, 'file: {}\nMean: {:.1f}\nstandard deviation: {:.1f}'.
              format(filename, lenmean, lensd), fontsize=10, transform=ax.transAxes)
    plt.axvline( lenmean, color='red', linewidth=1.5)

    plt.savefig('{}.pdf'.format(sample))
    plt.show()

    exit(0)

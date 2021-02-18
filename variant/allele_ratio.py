"""=================================================================================================
calculate depth and minor allele fraction along a sequence from mpileup file

Michael Gribskov     05 February 2021
================================================================================================="""
import math
import matplotlib.pyplot as plt
from mpileup import Mpileup

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    bins = 100

    mp = Mpileup()
    mp.open_file('../data/PrFi_10k.mpileup', 'r')

    position = []
    depthpos = []
    mafpos = []
    maf_histogram = [0 for _ in range(100)]

    freq = {}
    while mp.parse():
        mp.countchar()
        depth = mp.parsed['depth']
        depthplus = depth + 1.0
        for c in 'ACGTN*':
            freq[c] = (mp.count[c] + 0.167) / depthplus

        freq_order = sorted(freq.keys(), key=lambda x: freq[x], reverse=True)
        major = freq_order[0]

        if not depth or mp.count[major] == depth:
            continue

        minorfreq = 1.0 - freq[major]
        mafpos.append(minorfreq)
        mafh = math.floor(bins * minorfreq)
        maf_histogram[mafh] += 1

        position.append(mp.parsed['position'])
        depthpos.append(depth)
        print('{}\t{}\t{}\t{}\t{:0.3f}'.format(position[-1], major, depthpos[-1], mafh, mafpos[-1]))

    # plot
    plt.plot(position,mafpos)
    plt.show()

    plt.plot(range(0,100), maf_histogram)
    plt.show()


    exit(0)

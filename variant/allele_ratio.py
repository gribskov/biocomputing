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
    mp.open_file('C:/users/michael/Desktop/apple/out5.mpileup', 'r')

    position = []
    depthpos = []
    mafpos = []
    maf_histogram = [0 for _ in range(bins)]

    freq = {}
    n = 0
    next = 100
    while mp.parse():

        mp.countchar()
        depth = mp.parsed['depth']
        if depth < 20:
            continue

        n += 1
        # if n < next:
        #     continue
        #
        # next += 10000
        depthplus = depth + 1.0
        for c in 'ACGTN*':
            freq[c] = (mp.count[c] + 0.167) / depthplus

        freq_order = sorted(freq.keys(), key=lambda x: freq[x], reverse=True)
        major = freq_order[0]

            # continue

        minorfreq = 1.0 - freq[major]
        mafpos.append(minorfreq)
        mafh = math.floor(bins * minorfreq)
        # maf_histogram[mafh] += 1

        position.append(mp.parsed['position'])
        depthpos.append(depth)
        if minorfreq > 0.15:
            print('{}\t{}\t{}\t{}\t{:0.3f}'.format(position[-1], major, depthpos[-1], mafh, mafpos[-1]))
            maf_histogram[mafh] += 1
        # if n > 100000:
        #     break
        if not depth or mp.count[major] == depth:
            continue
    # plot
    plt.plot(position,mafpos)
    plt.show()

    plt.plot(range(0,bins), maf_histogram)
    plt.show()


    exit(0)

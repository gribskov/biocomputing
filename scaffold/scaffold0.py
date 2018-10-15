"""-------------------------------------------------------------------------------------------------
scaffold

15 October 2018     Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys
import pysam

# --------------------------------------------------------------------------------------------------\
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    bamfile = sys.argv[1]
    sys.stderr.write('bam file: {}\n'.format(bamfile))

    bam = pysam.AlignmentFile(bamfile, 'rb')

    # ref = {}
    # nref = 0
    # for i in range(len(bam.references)):
    #     ref[bam.references[i]] = bam.lengths[i]
    #     print( 'ref:{}     {}'.format(bam.references[i], bam.lengths[i]))
    #     nref += 1
    #     if nref > 99:
    #         break

    sys.stderr.write('references: {}\n'.format(bam.nreferences))

    nbam = 0
    nref = 0
    for ref in bam.references:
        nref += 1
        if nref > 4:
            break

        p = 0
        print(ref)
        for pos in bam.pileup(ref):
            p += 1
            print(pos.n)
            if p > 5:
                break

        nread = 0
        for read in bam.fetch(ref):
            # print(read)
            nread += 1
            # if nread > 5:
            #     break

        print('nread:',nread)

exit(0)

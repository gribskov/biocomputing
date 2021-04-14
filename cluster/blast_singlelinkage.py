"""=================================================================================================
Single linkage clustering for all-against-all blast search

Michael Gribskov     12 April 2021
================================================================================================="""
import sys
from blast import Blast
from cluster.single_linkage import SingleLinkage

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    blast = Blast(sys.argv[1])
    fmt = 'qname sname id alignlen mismatch gapopen qbeg qend sbeg send evalue bit_score'
    nfields = blast.setFormat(fmt)

    out = open('clusters80.out', 'w')

    cluster = SingleLinkage()
    cluster.set_keys(['sname', 'qname', 'evalue', 'id'])
    cluster.labels = ['sname', 'qname']

    nseq = 0
    while blast.readTabular():
        nseq += 1
        b = {'sname': blast.sname, 'qname':blast.qname,
             'evalue':float(blast.evalue), 'id':float(blast.id)}
        cluster.append(b)

        if not nseq % 100000:
            print('.', end='')
        if not nseq % 5000000:
            print()

    print('{} hits read'.format(nseq))

    cluster.single(threshold=1e-80)
    count = cluster.write_groups(out, names=True)

    out.write('{} hits clustered\n'.format(nseq))
    out.write('{} sequences\n'.format(len(cluster.names)))
    out.write('groups size histogram (number of sequences, number of groups\n')
    for c in range(len(count)):
        if count[c] > 0:
            out.write('{}\t{}\n'.format(c, count[c]))


    exit(0)

"""=================================================================================================
pacbio_sort produces a list of pacbio sequences sorted by length. Each line has statistics about
the sequences in sorted order

>m64165_221224_025950/61735974/5026_187590                  3   182565  592153  0.0001
id                                                       rank   length cum.len  frac
where frac cum.length/total.length

Michael Gribskov     10 July 2024
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    pb = open(sys.argv[1], 'r')
    order = open(sys.argv[2], 'r')
    cutfile = sys.argv[2] + '.cut'
    cut = open(cutfile, 'w')

    cutoff_len = 15000

    # read the statistics
    keep = []
    n_read = 0
    for line in order:
        n_read += 1
        id, rank, length, cumulative, fraction = line.rstrip().split('\t')
        if length > cutoff_len:
            keep.append( id )

    order.close()
    print(f"{n_read} sequence read from sys.argv[1]")
    print(f"{len(keep)} will be kept\n")

    select = False
    n_select = 0
    n_base = 0
    for line in pb:
        if line.startswith('>'):
            select = False
            if line in keep:
                select = True
                n_select += 1
                cut.write(line)
        elif select:
            cut.write(line)
            n_base += len(line) - 1

    pb.close()
    cut.close()
    print(f"{n_select} sequences with {n_base} bases written to {cutfile}")

    exit(0)

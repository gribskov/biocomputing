"""=================================================================================================
sort pacbio CLR split reads by length

header lines look like
>m64165_221224_025950/11/1637_22321
the last field gives the begin_end position in the read, hence the length

Michael Gribskov     10 July 2024
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    pacbiofile = sys.argv[1]
    sortedfile = pacbiofile + '.order'
    pb = open(pacbiofile, 'r')

    reads = []

    total_length = 0
    for line in pb:
        if line.startswith('>'):
            id = line.rstrip()
            fields = id.split('/')
            begin, end = fields[-1].split('_')
            length = int(end) - int(begin) + 1
            reads.append({'id': id, 'length': length})
            total_length += length


    out = open(sortedfile, 'w')
    sum_length = 0
    n = 0
    for r in sorted(reads, key=lambda l: l['length'], reverse=True):
        n += 1
        sum_length += r['length']
        fraction = sum_length / total_length
        print(f"{r['id']}\t{n:5d}\t{fraction:.4f}\t{r['length']}\t{sum_length}")

    exit(0)

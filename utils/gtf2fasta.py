#!/usr/bin/env python3
"""-------------------------------------------------------------------------------------------------
Read feature ranges from GTF file and extract fast formatted sequence

usage
    gtf2fasta <gtf_file> <fasta_file>
-------------------------------------------------------------------------------------------------"""
import sys
from sequence.fasta import Fasta

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # open files
    gtffile = sys.argv[1]
    try:
        gtf = open(gtffile, 'r')
    except:
        sys.stderr.write('Unable to open GTF file ({})\n'.format(gtffile))
        exit(1)

    fastafile = sys.argv[2]
    try:
        fasta = open(fastafile, 'r')
    except:
        sys.stderr.write('Unable to open fasta file ({})\n'.format(fastafile))
        exit(2)

    sys.stderr.write('gtf2fasta\n')
    sys.stderr.write('\tGTF: {}\n'.format(gtffile))
    sys.stderr.write('\tFasta: {}\n'.format(fastafile))

    features = ('gene', 'pseudogene' )
    nline = 0
    nfeature = 0
    for line in gtf:

        if line.startswith('#'):
            # skip comment lines
            continue

        print(line)
        nline += 1
        if nline > 50:
            break

        field = line.rstrip().split(maxsplit=8)
        info = {'seqname': field[0],
                'source': field[1],
                'feature': field[2],
                'begin': field[3],
                'end': field[4],
                'score': field[5],
                'strand': field[6],
                'frame': field[7],
                'attributes': field[8]}

        a = info['attributes'].split(';')
        attr = {}
        for note in a:
            # print('\tnote:{}:'.format(note))
            try:
                tag, value = note.split(maxsplit=1)
            except:
                continue
            # print('\ttag:{}:\tvalue:{}:'.format(tag, value))
            info[tag] = value.strip('"')
        del info['attributes']

        if info['feature'] in features:
            flist[info['ID']] = info
        for k in info:
            print('\t{}\t{}'.format(k, info[k]))

        print('')

exit(0)

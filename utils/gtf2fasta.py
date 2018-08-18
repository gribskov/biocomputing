#!/usr/bin/env python3
"""-------------------------------------------------------------------------------------------------
Read feature ranges from GTF file and extract fast formatted sequence

usage
    gtf2fasta <gtf_file> <fasta_file>
-------------------------------------------------------------------------------------------------"""
import sys
from sequence.fasta import Fasta


def complement(seq):
    """---------------------------------------------------------------------------------------------
    reverse complement the sequence
    ---------------------------------------------------------------------------------------------"""
    r = str.maketrans('ACGT', 'TGCA')
    return seq.translate(r)[::-1]


def formatseq(seq, linelen=100):
    """---------------------------------------------------------------------------------------------
    return a string with no more than linelen characters per line
    ---------------------------------------------------------------------------------------------"""
    fstr = ''
    begin = 0
    end = len(seq)
    while begin < end:
        fstr += seq[begin:begin + linelen] + '\n'
        begin += linelen

    return fstr


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

    seq = {}
    fasta = Fasta()
    fasta.open(sys.argv[2])
    sys.stderr.write('Reading Fasta {}...\n'.format(sys.argv[2]))
    nseq = 0
    while fasta.next():
        seq[fasta.id] = fasta.seq
        nseq += 1

    sys.stderr.write('\n{} Sequences read from {}\n'.format(nseq, sys.argv[2]))
    for s in seq:
        sys.stderr.write('\t{} len={}\n'.format(s, len(seq[s])))

    sys.stderr.write('\ngtf2fasta\n')
    sys.stderr.write('\tGTF: {}\n'.format(gtffile))
    sys.stderr.write('\tFasta: {}\n'.format(sys.argv[2]))

    features = ('gene', 'pseudogene', 'tRNA', 'rRNA')
    skip = ['region', 'direct_repeat', 'riboswitch']
    save = ['begin', 'end', 'strand', 'Name', 'Dbxref', 'gene_biotype', 'product', 'locus_tag']
    flist = {}
    nline = 0
    nfeature = 0
    for line in gtf:

        if line.startswith('#'):
            # skip comment lines
            continue

        sys.stderr.write('{}\n'.format(line))
        nline += 1
        # if nline > 50:
        #     break

        field = line.rstrip().split('\t', maxsplit=8)
        info = {'seqname': field[0],
                'source': field[1],
                'feature': field[2],
                'begin': int(field[3]),
                'end': int(field[4]),
                'score': field[5],
                'strand': field[6],
                'frame': field[7],
                'attributes': field[8]}

        # print('n:{} f:{}'.format(nline, info['feature']))
        if info['feature'] in skip:
            continue

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
            nfeature += 1
        elif info['Parent'] in flist:
            for k in info:
                if k not in flist[info['Parent']]:
                    flist[info['Parent']][k] = info[k]
        else:
            # flist[info['ID']] = info
            sys.stderr.write('unknown feature {}\n'.format(info['feature']))

     # write out sequences
    for gene in flist:
        thisgene = flist[gene]
        f = Fasta()
        f.id = thisgene['ID']
        f.doc = ''
        for k in save:
            if k in thisgene:
                f.doc += ' {}:{}'.format(k, thisgene[k])
        f.seq = seq[thisgene['seqname']][thisgene['begin'] - 1:thisgene['end']]
        if (thisgene['end'] - thisgene['begin'] > 100000):
            # coordinates cross origin
            f.seq = seq[thisgene['seqname']][thisgene['end'] - 1:] + seq[thisgene['seqname']][:thisgene['begin']]

        if thisgene['strand'] == '-':
            f.seq = complement(f.seq)

        sys.stdout.write(f.format(linelen=100))

exit(0)

"""=================================================================================================
Use the selected genes from DESeq2 prefiltering to select a subset of predicted transcripts from the
trinity output

write.table(rownames(filtered),file="filtered.names.txt", row.names=F, quote=F, col.names=F)

Michael Gribskov     12 November 2022
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    split_level = int(sys.argv[1])
    genefile = sys.argv[2]
    trinityfile = sys.argv[3]
    if len(sys.argv) > 4:
        outfile = sys.argv[4]
    else:
        outfile = 'selected.fa'
    out = open(outfile, 'w')

    # list of names, possibly at any level, cluster, component, gene, isoform
    sys.stderr.write(f'filtered gene file: {genefile}\n')
    gene = open(genefile, 'r')

    genes = []
    genecol = 1
    prefix = 'TRINITY_'
    for line in gene:
        if not line:
            continue

        if line.find('\t') > -1:
            name = line.split('\t')[genecol]
            name = name.replace(prefix, '')
        else:
            name = line.rstrip().replace(prefix, '')

        genes.append(name)

    gene.close()
    sys.stderr.write(f'genes found: {len(genes)}\n')

    # read through trinity fasta file, keeping all sequences that begin with a name in the list
    sys.stderr.write(f'trinity fasta file: {trinityfile}\n')
    trinity = open(trinityfile, 'r')

    name = ''
    rest = ''
    nout = 0
    found = []
    splitlevel = 5
    shortlevel = 3
    ntranscript = 0
    for line in trinity:
        if line.startswith('>'):
            # beginning of sequence, check name
            ntranscript += 1
            field = line.split('_')
            lfield = field[4].split(' ')
            del field[4]
            field += lfield
            name = '_'.join(field[1:splitlevel])
            keep = name in genes
            if keep:
                nout += 1
                if name not in found:
                    found.append(name)

                shortname = '_'.join(field[1:shortlevel])
                out.write(f">{shortname} {name} {field[5]}\n")
                # if not nout % 1000:
                #     print('.', end='')
                # if not nout % 100000:
                #     print(f' {nout}')
                sys.stderr.write(f'\r{nout}/{ntranscript}')
        else:
            if keep:
                out.write(line)

    sys.stderr.write(f'\ngenes found: {len(found)}\n')
    sys.stderr.write(f'transcripts written: {nout}\n')

    exit(0)

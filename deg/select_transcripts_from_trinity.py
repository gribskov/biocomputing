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
    # prefix = '>TRINITY_'
    for line in gene:
        if not line:
            continue
        genes.append(line.rstrip())

    gene.close()
    sys.stderr.write(f'genes found: {len(genes)}\n')

    # read through trinity fasta file, keeping all sequences that begin with a name in the list
    sys.stderr.write(f'trinity fasta file: {trinityfile}\n')
    trinity = open(trinityfile, 'r')

    name = ''
    rest = ''
    nout = 0
    found = []
    for line in trinity:
        if line.startswith('>'):
            # beginning of sequence, check name
            field = line.split('_')
            name = '_'.join(field[1:split_level])
            keep = name in genes
            if keep:
                nout += 1
                if name not in found:
                    found.append(name)
                shortname = ' '.join(line.split(' ')[:2])
                shortname = shortname.replace('>TRINITY_','')
                out.write(f">{shortname}\n")
                # if not nout % 1000:
                #     print('.', end='')
                # if not nout % 100000:
                #     print(f' {nout}')
                sys.stderr.write(f'\r{nout}')
        else:
            if keep:
                out.write(line)

    sys.stderr.write(f'\ngenes found: {len(found)}\n')
    sys.stderr.write(f'transcripts written: {nout}\n')

    exit(0)

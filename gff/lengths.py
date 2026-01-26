"""#################################################################################################
lengths.py
get the intergenic region and intron lengths based on gff file
use gene and exon features, respectivelyh

2025-11-25 gribskov
#################################################################################################"""
from gff2 import GxfRecord
from gff2 import GxfSet

####################################################################################################
# Main
####################################################################################################
if __name__ == '__main__':
    gtffile = 'data/Plantago_major_annotation.gff'
    features = ['gene', 'exon']

    print(f'lengths.py')
    print(f'\tGTF file: {gtffile}')
    print(f'\tSelected features: {features}')

    # read in gtf with selected features
    gtf = GxfSet(file=gtffile, fmt='gff')
    feature_n = gtf.feature_get(features)
    print(f'\nfeatures read from {gtffile}: {feature_n} ')

    glen = []
    elen = []
    seq = ''

    for f in sorted(gtf.features, key=lambda x: (x.seqid, x.start)):
        if f.seqid != seq:
            gene_old = ''
            gene_end = 0
            exon_old = ''
            exon_end = None
            seq = f.seqid

        if f.type == 'gene':
            gene = f.attribute['ID']
            gene = gene.replace('-gene', '')

            gl = 0
            if gene_old:
                gl = f.start - gene_end
                glen.append(gl)
            gene_old = gene
            gene_end = f.end

            exon = ''
            exon_old = ''
            exon_end = None
            # print(f'gene: {gene}\t{f.start}\t{f.end}\t{f.strand}\t{gl}')

        if f.type == 'exon':
            el = 0
            if exon_old:
                el = f.start - exon_end
                elen.append(el)
            exon_old = f.attribute['Parent']
            exon_end = f.end
            # print(f'exon: {exon_old}\t{f.start}\t{f.end}\t{f.strand}\t{el}')


    for l in elen:
        print(f'{l}')

    exit(0)

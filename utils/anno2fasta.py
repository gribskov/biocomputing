####################################################################################################
# written for hass avocado v2 genome annotation
# data starts at row 6
# files are tab separated with fields
# 0 blank column
# 1 row
# 2 Gene ID
# 3 Length (nt)
# 4 cDNA Sequences (nt)
# 5 Length(aa)
# 6 Protein Sequences (aa)
# 7% of the homologous proteins represented by avocado gene models (Average)
# 8 Status for the Predicted Gene Model Partial/Complete
# 9 blank column
# 10-18 Annotation based on top-BLAST-hit method Amborella trichopoda K-S
#     Accession No.
#     Description
#     Length(aa)
#     Protein Sequence
#     Identity %
#     Alignment length
#     e-value
#     Score
#     blank column
# 19-27 Annotation based on top-BLAST-hit method Sorghum bicolor T-AA
#     Accession No.
#     Description
#     Length(aa)
#     Protein Sequence
#     Identity %
#     Alignment length
#     e-value
#     Score
#     blank column
# 28-36 Annotation based on top-BLAST-hit method Vitis vinifera AC-AJ
#     Accession No.
#     Description
#     Length(aa)
#     Protein Sequence
#     Identity %
#     Alignment length
#     e-value
#     Score
#     blank column
# 37-45 Annotation based on top-BLAST-hit method Solanum lycopersicum AL-AT
#     Accession No.
#     Description
#     Length(aa)
#     Protein Sequence
#     Identity %
#     Alignment length
#     e-value
#     Score
#     blank column
#  46-55 Annotation based on top-BLAST-hit method Arabidopsis thaliana AU-BD
#     Accession No.
#     Symbol
#     Description
#     Length(aa)
#     Protein Sequence
#     Identity %
#     Alignment length
#     e-value
#     Score
#     blank column
#
# 56-61 Pfam domains BE-BJ
#     number of Pfam domain identified
#     e-value
#     score
#     Accesion
#     Name
#     blank column
#
# 62-66 Funtional annotation based on KEGG BK-BO
#     KO terms
#     Symbol
#     Name
#     Enzyme Commission (EC) number
#     blank column
#
# 67-73 Funtional annotation based on GO BP-BU
#     number of cellular components
#     cellular components GO terms
#     number of molecular functions
#     molecular functions GO terms
#     number of biological processes
#     biological processes GO terms
####################################################################################################
import sys

# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    infile = sys.argv[1]
    sys.stderr.write('Blast search: {}\n'.format(infile))
    try:
        anno = open(infile, 'r')
    except OSError:
        sys.stderr.write('Cannot open annotationfile ({})'.format(infile))
        exit(1)

    # output files
    pfile = '{}.pfa'.format(infile)
    sys.stderr.write('\tProtein output: {}\n'.format(pfile))
    try:
        pro = open(pfile, 'w')
    except OSError:
        sys.stderr.write('Cannot open protein output file ({})\n'.format(pfile))
        exit(1)

    nfile = '{}.fa'.format(infile)
    sys.stderr.write('\tcDNA output: {}\n'.format(nfile))
    try:
        nuc = open(nfile, 'w')
    except OSError:
        sys.stderr.write('Cannot open cDNA output file ({})\n'.format(nfile))
        exit(1)

    # skip header
    nline = 0
    for line in anno:
        # print('{} {}'.format(nline,line))
        nline += 1
        if nline > 5:
            break

    # main data
    nline = 0
    linelen = 60
    for line in anno:
        if not line:
            continue

        field = line.rstrip().split('\t')
        try:
            # there are unparseable notes at the end of the file
            nuclen = int(field[2])
        except IndexError:
            continue

        id = field[1].replace('|HassGeneModel','')
        nucseq = field[3]
        nuc.write('>{} len={}\n'.format(id, nuclen))
        pos = 0
        while pos < nuclen:
            nuc.write('{}\n'.format(nucseq[pos:pos + linelen]))
            pos += linelen

        prolen = int(field[4])
        proseq = field[5]
        pro.write('>{} len={}\n'.format(id, prolen))
        pos = 0
        while pos < prolen:
            pro.write('{}\n'.format(proseq[pos:pos + linelen]))
            pos += linelen

        nline += 1
        # if nline > 5:
        #     break

    sys.stderr.write('\t{} gene models written\n'.format(nline))

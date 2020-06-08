"""=================================================================================================
merge kraken2 reports into a table and filter by number of samples, counts of taxa present

 64.95  1895470 1895470 U       0       unclassified
 35.05  1022732 874     R       1       root
 34.34  1001967 29770   R1      131567    cellular organisms
 22.72  663149  34      D       2157        Archaea
 22.71  662861  95      P       28890         Euryarchaeota
 22.70  662418  236     P1      2290931         Stenosarchaea group
 22.67  661544  51822   C       183963            Halobacteria
 10.70  312220  13140   O       2235                Halobacteriales
  6.19  180549  9776    F       1963268               Haloarculaceae
  1.61  47094   1278    G       63743                   Natronomonas
  0.85  24773   0       S       416273                    Natronomonas moolapensis
  0.85  24773   24773   S1      268739                      Natronomonas moolapensis 8.8.11


columns:
1. % of minimizers mapping to clade
2. number of minimizers mapping to clade
3. number of minimizers mapping to taxon
4. taxonomic rank, note that integers may be appended to indicate subranks
    U unclassified
    R root
    D domain (kingdom)
    P phylum
    C class
    O order
    F family
    G genus
    S species
5. taxon ID
6. common name of taxon (indentation is significant)


Michael Gribskov     07 June 2020
================================================================================================="""
import sys

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # open file
    # idlist
    reportname = sys.argv[1]
    sys.stderr.write('\tKraken2 report: {}\n'.format(reportname))
    try:
        report = open(reportname, 'r')
    except:
        sys.stderr.write('\nUnable to open Kraken2 report ({})\n'.format(reportname))
        exit(1)

    # define rank order for taxonomic ranks
    i2r = ['U', 'D', 'P', 'C', 'O', 'F', 'G', 'S' ]
    for i in range(len(i2r)):
        r2i[i2r[i]] = i

    # read and store in array of hashes
    taxonomy = []
    for line in report:
        ( pct_mapped, n_mapped, n_taxon, rank, taxid, text ) = line.rstrip().split('\t')
        print('{}: {}'.format(pct_mapped, text))
        taxonomy.append(    {   'pct_mapped':pct_mapped,
                                'n_mapped':n_mapped,
                                'n_taxon':n_taxon,
                                'rank':rank,
                                'text':text
        })

    exit(0)
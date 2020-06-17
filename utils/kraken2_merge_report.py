"""=================================================================================================
merge kraken2 reports into a table both raw and normalized counts are shown

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
import os.path
import glob

from taxonomy.taxonomy import *


def openReport(reportname):
    """---------------------------------------------------------------------------------------------
    Safe open report file for reading, return filehandle.  Unopenable files are skipped with
    error message

    :param reportname: string, readable file name
    :return: open filehandle
    ---------------------------------------------------------------------------------------------"""
    try:
        report = open(reportname, 'r')
    except:
        sys.stderr.write('\nopenReport - Unable to open Kraken2 report ({})\n'.format(reportname))

    return report


def writeall(file, mergedtax, taxset, sep=''):
    """---------------------------------------------------------------------------------------------
    Write the final result 
    
    :param file: filehandle for output
    :param mergedtax: Taxonomy with merged result
    :param taxset: list of Taxonomy with sample results
    :return: True
    ---------------------------------------------------------------------------------------------"""
    file.write('#\t{:>8s}\t{:>6s}\t{:>6s}'.format('', 'Merged', ''))
    for sample in taxset:
        file.write('{}{:>8s}{:>6s}\t{:>6s}'.format(sep, '', sample, ''))
    file.write('\n')

    file.write('#\t{:>8s}\t{:>6s}\t{:>6s}'.format('pct', 'clade', 'taxon'))
    for sample in taxset:
        file.write('{}{:>8s}\t{:>6s}\t{:>6s}'.format(sep, 'pct', 'clade', 'taxon'))
    file.write('{}{:>4}\t{:>6s}\t{}\n'. \
               format(sep, 'taxa', 'taxid', 'classification'))

    for node in mergedtax:
        taxid = node.taxid
        thisnode = mergedtax.index[taxid]
        file.write('\t{:8.2f}\t{:6d}\t{:6d}'. \
                   format(thisnode.pct_mapped, int(thisnode.n_mapped), int(thisnode.n_taxon)))

        ntaxa = 0
        for sample in taxset:
            try:
                thisnode = taxset[sample].index[taxid]
                pct_mapped = thisnode.pct_mapped
                ntaxa += 1
            except KeyError:
                # this taxonomy does not have this taxid so its count is zero
                pct_mapped = 0.0
            file.write('{}{:8.2f}\t{:6d}\t{:6d}'. \
                       format(sep, pct_mapped, int(thisnode.n_mapped), int(thisnode.n_taxon)))

        level = Taxonomy.r2i[node.rank[0]]
        space = ' ' * level
        file.write('{}{:>4d}\t{:>6d}\t{}{}\n'. \
                   format(sep, ntaxa, taxid, space, mergedtax.index[taxid].text))

    return True


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    # first argument is a multiple file spec
    target = sys.argv[1]

    # store taxonomy results in an array of dicts indexed by the file name
    taxset = {}

    nreport = 0
    sys.stderr.write('Kraken2 reports matching {}\n\n'.format(target))
    for reportname in glob.glob(target):
        nreport += 1

        report = openReport(reportname)
        sample = os.path.basename(reportname)
        find = sample.rfind('.')
        if find > -1:
            sample = sample[:find]
        sys.stderr.write('\t{} {} ... '.format(nreport, sample))

        # read and store in taxonomy object
        taxset[sample] = Taxonomy.readKraken(report)
        sys.stderr.write('{} taxa read\n'.format(len(taxset[sample].index)))

    # create merged taxonomy
    sys.stderr.write('\nMerging {} taxonmies\n\n'.format(nreport))
    mergedtax = Taxonomy()
    for tax in taxset:
        sys.stderr.write('\tMerging {} ... '.format(tax))
        ntaxa = mergedtax.merge(taxset[tax])
        sys.stderr.write('{} taxa after merging\n'.format(ntaxa))

    sys.stderr.write('\nMerged taxonmy\n\n')
    fmt = mergedtax.formatString()
    mergedtax.dumpFromTree(file=sys.stderr, fmt=fmt)

    # Write out merged result
    totalread = mergedtax.recalcPercent()
    sys.stdout.write('# total reads analyzed:{}\n'.format(totalread))
    sys.stdout.write('# Unnormalized counts\n')
    writeall(sys.stdout, mergedtax, taxset)

    # normalized count, rescale to integer count with n_mapped = target for classified reads
    # the sum of classified reads is in tax.root.child[0]

    target = 100000

    # two scale factors:
    #   scale: corrects the classified counts to equal target (integer rounded)
    #   ratio: corrects percentages to be percent of classified
    scale = {}
    ratio = {}
    nsample = len(taxset)
    for sample in taxset:
        unclassified = taxset[sample].root
        classified = taxset[sample].root.child[0]
        scale[sample] = target / classified.n_mapped
        ratio[sample] = (unclassified.n_taxon + classified.n_mapped) / classified.n_mapped

    for node in mergedtax:
        taxid = node.taxid
        merged = mergedtax.index[taxid]
        merged.pct_mapped = 0.0
        merged.n_mapped = 0
        merged.n_taxon = 0
        merged.ntaxa = 0

        for sample in taxset:
            try:
                thisnode = taxset[sample].index[taxid]
                thisnode.n_taxon = int(scale[sample] * thisnode.n_taxon)
                thisnode.n_mapped = int(scale[sample] * thisnode.n_mapped)
                thisnode.pct_mapped = ratio[sample] * thisnode.pct_mapped
                merged.pct_mapped += thisnode.pct_mapped
                merged.n_mapped += thisnode.n_mapped
                merged.n_taxon += thisnode.n_taxon


            except KeyError:
                pass

        merged.pct_mapped /= nsample
        merged.n_mapped /= nsample
        merged.n_taxon /= nsample

    sys.stdout.write('\n# Normalized counts: target={}\n'.format(target))
    writeall(sys.stdout, mergedtax, taxset, sep='  |')

    exit(0)

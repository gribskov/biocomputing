"""=================================================================================================
Given merged salmon counts (rows: trinity assemblies, columns: samples) merge at a specific
trinity level (bundle, component, gene, isoform) and select genes that are taxonomically good
(list produced by blast/blast_to_taxonomy) and that have at least a minimum number of counted reads
after merging.

Michael Gribskov     18 March 2025
================================================================================================="""


class Cluster:
    """---------------------------------------------------------------------------------------------
    Cluster is a group of Trinity assemblies that have been defined to be "the same"
    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        each sequence in seqlist is {'id', 'length', 'evalue', 'tax', 'lineage', 'good' }
        -----------------------------------------------------------------------------------------"""
        self.seqlist = []
        self.ngood = 0
        self.nbad = 0


def whitelist_get(goodfile):
    """---------------------------------------------------------------------------------------------
    Read the file written by blast/blast_to_taxonomy format is
    51: rosids	Eukaryota;Viridiplantae;Streptophyta;Magnoliopsida
    skip lines that start with !

    :param goodfile: string     path to good taxonomy file
    :return: dict
    ---------------------------------------------------------------------------------------------"""
    good = open(goodfile, 'r')
    good_taxa = {}
    for line in good:
        if line.startswith('!') or not line.rstrip():
            continue

        n, tax, lineage = line.rstrip().split(maxsplit=2)
        good_taxa[tax] = lineage

    good.close()
    return good_taxa


def blast_read(blastfile):
    """---------------------------------------------------------------------------------------------
    Read the blast result and store as a list of cluster objects

    :param blastfile:
    :return:
    ---------------------------------------------------------------------------------------------"""
    column = ['qid', 'qlen', 'qbegin', 'qend',
              'sid', 'slen', 'sbegin', 'send',
              'allen', 'pident', 'score', 'evalue', 'stitle']
    column_n = len(column)

    blast = open(blastfile, 'r')
    cluster = {}
    qid_old = ''
    blast = open(blastfile, 'r')
    nline = 0
    for line in blast:
        nline += 1

        field = line.rstrip().split('\t', maxsplit=column_n - 1)
        hit = {column[i]: field[i] for i in range(len(field))}
        print(f"{nline:7d}\t{hit['qid']}\t{hit['sid']}")

        qid = hit['qid']
        if qid != qid_old:
            cluster[qid] = Cluster()
            this = cluster[qid]
            qid_old = qid

        taxname = blast_get_tax(hit['stitle'])
        this.seqlist.append({'length':  int(hit['qlen']),
                             'evalue':  float(hit['evalue']),
                             'tax':     taxname,
                             'lineage': '',
                             'good':    0
                             })

    blast.close()
    return cluster


def blast_get_tax(stitle):
    """---------------------------------------------------------------------------------------------
    get the Tax string from the subject title field

    :param stitle: string       subject title from blast search
    :return: string             Tax string
    ---------------------------------------------------------------------------------------------"""
    taxpos = stitle.find('Tax=') + 4
    taxidpos = stitle.find('TaxID')
    tax = stitle[taxpos:taxidpos]

    # filter to match the output of blast_to_taxonomy.py
    if tax.find('(') != -1:
        # trim off information in parentheses
        # this information is usually not part of the actual taxonomic name
        tax = tax[:tax.find('(') - 1]

    # remove strain information (sometimes causes multiple returns from taxonomy server)
    strain = tax.find('_str.')
    if strain > -1:
        tax = tax[:strain]

    # remove terms uncultured and unclassified
    uncultured = tax.find('uncultured_')
    if uncultured > -1:
        tax = tax[uncultured + 11:]

    unclassifed = tax.find('unclassified_')
    if unclassifed > -1:
        tax = tax[unclassifed + 13:]

    # convert spaces and / to underline
    tax = tax.strip().replace(' ', '_')
    tax = tax.replace('/', '_')

    return tax


def add_tax_info(trinity, goodtax):
    """---------------------------------------------------------------------------------------------
    add the taxonomy information from goodtax to the list of trinity clusters. Note that only "good"
    assemplies will have lineage information

    :param trinity: list of Cluster     The blast information for each trinity assembly
    :param goodtax: dict                keys should correspond to tax in the trinity results
                                        values are the lineage
    :return: int                        number added
    ---------------------------------------------------------------------------------------------"""
    ngood = 0
    for iso in trinity:
        for seq in trinity[iso].seqlist:
            tax = seq['tax']
            if tax in goodtax:
                seq['lineage'] = goodtax[tax]
                seq['good'] = 1
                ngood += 1
    return ngood

def merge_at_level(trinity, level):
    """---------------------------------------------------------------------------------------------
    merge the isoform data at the selected level. 2=bundle, 3=component, 4=gene

    :param trinity: list of Cluster     blast and lineage data for each isoform
    :param level: int                   merging level (see level dict in main)
    :return: list of Cluster            clustered data
    ---------------------------------------------------------------------------------------------"""
    merged = {}
    for iso in trinity:
        # get the clustered ID
        field = iso.split('_')
        mid = '_'.join(field[:level])
        if mid not in merged:
            merged[mid] = Cluster()
            this = merged[mid]
        this.seqlist += trinity[iso].seqlist
        for seq in trinity[iso].seqlist:
            if seq['good']:
                this.ngood += 1
            else:
                this.nbad += 1

    return merged


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    goodfile = '../blast/data/c16c31.trinity.goodtax.txt'
    blastfile = '../blast/data/c16c31.trinity_uniref_1e-10.dmndblastx'
    countfile = 'data/C16C31.numreads.txt'

    # isoform would be no merging
    level = {'bundle':2, 'component':3, 'gene':4}

    # read in taxonomy white list
    good_taxa = whitelist_get(goodfile)
    print(f'Good taxa: {len(good_taxa)} read from {goodfile}')

    # read blast search and characterize each query (trinity assembly) as good or bad based on
    # multiple hits
    trinity = blast_read(blastfile)
    print(f'trinity isoforms: {len(trinity)} read from {blastfile}')

    nannot = add_tax_info(trinity, good_taxa)
    print(f'transcripts annotated: {nannot}')
    merged = merge_at_level(trinity, level['gene'])
    print(f'merged assemblies at gene level: {len(merged)}')

    for gene in merged:
        this = merged[gene]
        if merged[gene].ngood and merged[gene].nbad:
            # ambiguous?
            print('ambiguous')

    # write out with total count cutoff
    prefix = 'counts.merged.'
    # open files for each output
    total = []
    count = open(countfile, 'r')
    header = count.readline()
    for line in count:


    exit(0)

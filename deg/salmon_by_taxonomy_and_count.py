"""=================================================================================================
Given merged salmon counts (rows: trinity assemblies, columns: samples) merge at a specific
trinity level (bundle, component, gene, isoform) and select genes from labelled taxonomic groups
(list produced by blast/blast_to_taxonomy) and that have at least a minimum number of counted reads
after merging.

Michael Gribskov     18 March 2025
================================================================================================="""
from collections import defaultdict


class Cluster:
    """---------------------------------------------------------------------------------------------
    Cluster is a group of Trinity assemblies that have been defined to be "the same"
    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        each sequence in seqlist is {'id', 'length', 'evalue', 'tax', 'lineage', 'good' }
        -----------------------------------------------------------------------------------------"""
        self.seqlist = []
        self.ngroup = defaultdict(int)
        self.group = ''
        self.colcount = []
        self.total = 0
        total = 0

    def count_add(self, counts):
        """----------------------------------------------------------------------------------------
        Add the counts for each datacolumn for all constituent isoforms

        :param counts:
        :return:
        ----------------------------------------------------------------------------------------"""
        colcount = self.colcount
        countfloat = [float(count) for count in counts]
        if not colcount:
            # colcount is empty
            colcount = countfloat
        else:
            # add to existing counts
            for i in range(len(countfloat)):
                colcount[i] += countfloat[i]

        total = 0
        for i in range(len(countfloat)):
            total += countfloat[i]

        self.total += total
        self.colcount = colcount

        return total


def whitelist_get(goodfile):
    """---------------------------------------------------------------------------------------------
    No longer used
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


def groups_read(groupfile):
    """---------------------------------------------------------------------------------------------
    Read tax and groups from output of blast/blast_to_taxonomy.py
    Lines beginning with ! are skipped

    format
    1179       plant     Wenchengia_alternifolia	Eukaryota;Viridiplantae;Streptophyta;...

    :param groupfile: string    path to file with taxa/group information
    :return: dict               keys: taxon, values: dict of {'group', 'lineage'}
    ---------------------------------------------------------------------------------------------"""
    info = open(groupfile, 'r')
    groups = {}
    for line in info:
        if line.startswith('!') or not line.rstrip():
            continue

        try:
            n, group, taxon, lineage = line.rstrip().split()
        except ValueError:
            print(line)

        groups[taxon] = {'group': group, 'lineage': lineage}

    info.close()
    return groups


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
                             'group':   ''
                             })
        # if nline > 10000:
        #     # TODO remove debugging code
        #     break

    blast.close()
    return cluster


def blast_get_tax(stitle):
    """---------------------------------------------------------------------------------------------
    get the Tax string from the subject title field. This has to duplicate the filtering in
    blast/blast_to_taxonomy.py or the taxa won't match

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
    strain = tax.find(' str.')
    if strain > -1:
        tax = tax[:strain]

    # remove terms uncultured and unclassified
    uncultured = tax.find('uncultured ')
    if uncultured > -1:
        tax = tax[uncultured + 11:]

    unclassifed = tax.find('unclassified ')
    if unclassifed > -1:
        tax = tax[unclassifed + 13:]

    # convert spaces and / to underline
    tax = tax.strip().replace(' ', '_')
    tax = tax.replace('/', '_')

    return tax


def add_tax_info(trinity, groups):
    """---------------------------------------------------------------------------------------------
    add the taxonomy information from groups to the list of trinity clusters.

    :param trinity: list of Cluster     The blast information for each trinity assembly
    :param groups: dict                 keys should correspond to tax in the trinity results
                                        values are the lineage
    :return: dict                       keys are groups, values are counts of annotated entries
    ---------------------------------------------------------------------------------------------"""
    count = defaultdict(int)
    for iso in trinity:
        for seq in trinity[iso].seqlist:
            tax = seq['tax']
            seq['lineage'] = groups[tax]['lineage']
            seq['group'] = groups[tax]['group']
            count[seq['group']] += 1

    return count


def merged_id(trinity_id, level):
    """---------------------------------------------------------------------------------------------
    truncate the trinity at the selected level. levels: 2=bundle, 3=component, 4=gene

    :param trinity_id: string       trinity ID such as TRINITY_DN203973_c0_g1_i1
    :param level: int               merging level (see level dict in main)
    :return: string                 truncated ID such as TRINITY_DN203973_c0_g1
    ---------------------------------------------------------------------------------------------"""
    field = trinity_id.split('_')
    mid = '_'.join(field[:level])

    return mid


def merge_at_level(trinity, level):
    """---------------------------------------------------------------------------------------------
    merge the isoform data at the selected level. 2=bundle, 3=component, 4=gene

    :param trinity: list of Cluster     blast and lineage data for each isoform
    :param level: int                   merging level (see level dict in main)
    :return: list of Cluster            clustered data
    ---------------------------------------------------------------------------------------------"""
    merged = {}
    for iso in trinity:
        mid = merged_id(iso, level)

        if mid not in merged:
            merged[mid] = Cluster()
            this = merged[mid]

        this.seqlist += trinity[iso].seqlist
        for seq in trinity[iso].seqlist:
            merged[mid].ngroup[seq['group']] += 1

    return merged


def group_assign(merged):
    """---------------------------------------------------------------------------------------------
    assign the merged isoforms to a group based on the counts of the groups of each constituent
    isoform. The rules are:
    1) assign as plant if there are any plants
    2) assign as highest count
    3) if it's a tie, assign as other

    :param merged: list of Cluster      annotated merged clusters
    :return: string                     count of number assigned to each group
    ---------------------------------------------------------------------------------------------"""
    count = defaultdict(int)
    group = 'other'
    for m in merged:
        this = merged[m]
        if 'plant' in this.ngroup:
            this.group = 'plant'
            count['plant'] += 1
            continue

        gsort = sorted(this.ngroup, key=lambda g: this.ngroup[g], reverse=True)
        if len(gsort) == 1:
            gbest = gsort[0]
        elif gsort[0] > gsort[1]:
            gbest = gsort[0]
        else:
            gbest = 'other'

        this.group = gbest
        count[gbest] += 1

    return count


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    groupfile = '../blast/data/c16c31.trinity.goodtax.txt'
    blastfile = 'data/c16c31.trinity_uniref_1e-4.dmndblastx'
    countfile = 'data/C16C31.numreads.txt'

    # isoform would be no merging
    level = {'bundle': 2, 'component': 3, 'gene': 4}

    # read in taxonomy and groups
    print(f'Taxa read from {groupfile}')
    groups = groups_read(groupfile)
    gcount = defaultdict(int)
    for taxon in groups:
        group = groups[taxon]['group']
        gcount[group] += 1

    for g in sorted(gcount, key=lambda g: gcount[g], reverse=True):
        print(f'\t{g:10s}\t{gcount[g]}')

    # read blast search and characterize each query (trinity assembly) as good or bad based on
    # multiple hits
    trinity = blast_read(blastfile)
    print(f'\ntrinity isoforms: {len(trinity)} read from {blastfile}')

    nannot = add_tax_info(trinity, groups)
    print(f'\nblast hits used for annotation:')
    for g in nannot:
        print(f'\t{g:10s}\t{nannot[g]}')

    merged = merge_at_level(trinity, level['gene'])
    gcount = group_assign(merged)
    print(f'\nmerged assemblies at gene level: {len(merged)}')
    # for g in gcount:
    #     print(f'\t{g:10s}\t{gcount[g]}')
    gcount['unknown'] = 0  # to count isoforms that had no blast hits

    total = {g: 0 for g in gcount}
    total['all'] = 0  # overall total count across all samples and groups
    count = open(countfile, 'r')
    header = count.readline()

    for line in count:
        field = line.split()
        mid = merged_id(field[0], level['gene'])

        if mid not in merged:
            merged[mid] = Cluster()
            merged[mid].group = 'unknown'

        group = merged[mid].group
        count = merged[mid].count_add(field[1:])
        try:
            total[group] += count
        except KeyError:
            m = merged[mid]
            print('undef')

        total['all'] += count

    # write out with total count cutoff
    prefix = 'counts.merged.'
    # open files for each output
    gfile = {}
    for g in gcount:
        fname = prefix + g + '.count'
        gfile[g] = open(fname, 'w')

    hfield = header.split()
    newheader = [hfield[0]] + ['total'] + hfield[1:]
    # print(newheader)
    header = '\t'.join([hfield[0]] + ['total'] + hfield[1:])
    for f in gfile:
        gfile[f].write(f'{header}\n')

    print(f'\nRead counts by group ({countfile})')
    for g in total:
        print(f'\t{g:10s}\t{total[g]:>12.1f}')

    # write out counts to each group file
    n_written = 0
    print('\nwriting counts')
    for mid in sorted(merged, key=lambda c: merged[c].total, reverse=True):
        count = merged[mid].colcount
        group = merged[mid].group
        total = merged[mid].total
        fh = gfile[group]
        fh.write(f'{mid:30s}\t{total:11.1f}')
        for c in count:
            fh.write(f'\t{c:11.1f}')
        fh.write('\n')

        n_written += 1
        if not n_written % 1000:
            print('.', end='')
        if not n_written % 100000:
            print(f' {n_written}')
    print(f' {n_written}')
    print(f'\nCounts written to')
    for g in gfile:
        print(f'\t{gfile[g].name}')

    exit(0)

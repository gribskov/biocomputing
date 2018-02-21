"""=================================================================================================
Calculate the coverage for each trinity predicted transcript and tabulate by gene, component,
================================================================================================="""
import re
import sys
from blast.blast import Blast

idre = re.compile(r'>*TRINITY_([^_]+)_c(\d+)_g(\d+)_i(\d+)')  # regex for splitting trinity ID


def getIDParts(idstring):
    """---------------------------------------------------------------------------------------------
    Breakdown the trinity ID string to give the separate parts of the ID:
    Cluster,  component, gene and isoform
    :param:idstring trinity ID string to parse
    :return: dict of cluster,  component, gene, isoform, shortid
     --------------------------------------------------------------------------------------------"""
    if 'TRINITY' in idstring:
        id = {}
        cluster, component, gene, isoform = idre.match(idstring).groups()
        id['cluster'] = cluster
        id['component'] = component
        id['gene'] = gene
        id['isoform'] = isoform
        id['shortid'] = '{cl}.{co}.{g}.{i}'.format(cl=cluster, co=component, g=gene, i=isoform)

    else:
        # old trinity id string
        # comp1000022_c0_seq1
        component, gene, isoform = idstring.split('_')
        id['cluster'] = None
        id['component'] = component
        id['gene'] = gene
        id['isoform'] = isoform
        id['shortid'] = '{co}.{g}.{i}'.format(co=component, g=gene, i=isoform)

    return id

    # end of getIDParts


def histogram(coverage):
    """---------------------------------------------------------------------------------------------
    print a list with a histogram of the values at 1% intervals
    :param coverage: dictionary of coverage values
    :return: None
    ---------------------------------------------------------------------------------------------"""
    hist = [0 for i in range(101)]
    n = 0
    for key in sorted(coverage.values()):
        hist[int(key * 100)] += 1
        n += 1

    cum = 0
    for i in range(101):
        frac = hist[i] / n
        cum += frac
        print('\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}'.format(i, hist[i], frac, cum, 1.0 - cum))

    return None

    # end of histogram


# ==================================================================================================
# main program
# ==================================================================================================

# TODO: make this command line option
OLDID = True

# get input file name from command line
blastfilename = ''
if len(sys.argv) > 1:
    blastfilename = sys.argv[1]
else:
    print('Error: Please provide a BlastX/Diamond search result as input')
    exit(1)
print('  BlastX file:', blastfilename)

blast = Blast()
# blast.setFormat(
#     'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore')
blast.setFormat('qseqid qlen qstart qend sseqid slen sstart send length pident evalue stitle')

if not blast.new(blastfilename):
    print('Error opening blastX input file ({})'.format(blastfilename))
    exit(2)

nhits = 0
levels = ('cluster', 'component', 'gene', 'isoform')
if OLDID:
    levels = ('component', 'gene', 'isoform')

count = {k: {} for k in levels}
best = {k: {} for k in levels}
n = {k: 0 for k in levels}

# count the hits at each level
while blast.next():
    nhits += 1
    id = getIDParts(blast.qseqid)
    subj_cov = (int(blast.send) - int(blast.sstart) + 1) / int(blast.slen)
    item = ''
    for k in levels:
        item += id[k]
        try:
            if subj_cov > count[k][item]:
                count[k][item] = subj_cov
                best[k][item] = '_'.join(['%s' % id[k] for k in levels])
        except KeyError:
            # if item is undefined, save without testing
            count[k][item] = subj_cov
            best[k][item] = '_'.join(['%s' % id[k] for k in levels])
            n[k] += 1
        item += '_'

    # if nhits >= 100000:
    #     break

# print summary
print('{} blast hists processed'.format(nhits))
for k in levels:
    print('\t{} {}s'.format(n[k], k))

# generate and print histograms for each level
for k in levels:
    print('\n{} {}s'.format(n[k], k))
    histogram(count[k])

for k in levels:
    for i in count[k]:
        print('{}\t{}\t{:.3f}'.format(i, best[k][i], count[k][i]))

# print out best predicted transcripts
for k in levels:
    outfile = 'best_' + k + '.list'
    try:
        out = open(outfile, 'w')
    except:
        print('\nCould not open {} for writing'.format(outfile))
        continue

    print('Writing best {}'.format(k))
    for item in best[k]:
        out.write('{}\t{}\n'.format(item, best[k][item]))

    out.close()

exit(0)

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
    id = {}
    cluster, component, gene, isoform = idre.match(idstring).groups()
    id['cluster'] = cluster
    id['component'] = component
    id['gene'] = gene
    id['isoform'] = isoform
    id['shortid'] = '{cl}.{co}.{g}.{i}'.format(cl=cluster, co=component, g=gene, i=isoform)

    return id

    # end of getIDParts


# ==================================================================================================
# main program
# ==================================================================================================

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

cluster = {}
component = {}
gene = {}
isoform = {}
id_old = dict(cluster='', component='', gene='', isoform='', shortid='')
nhits = 0
ncluster = 0
ncomponent = 0
ngene = 0
nisoform = 0
while blast.next():
    nhits += 1
    # print('query:', blast.qseqid, 'subject:', blast.sseqid)
    id = getIDParts(blast.qseqid)
    subj_cov = (int(blast.send) - int(blast.sstart) + 1) / int(blast.slen)
    # print('  {:.3f}'.format(subj_cov))
    # print('  {cluster}  {component}  {gene}  {isoform}  short_id:{shortid}'.format(**id))

    key = id['cluster']
    try:
        if subj_cov > cluster[key]:
            cluster[key] = subj_cov
    except KeyError:
        cluster[key] = subj_cov
        ncluster += 1

    key += '_{}'.format(id['component'])
    try:
        if subj_cov > component[key]:
            component[key] = subj_cov
    except KeyError:
        component[key] = subj_cov
        ncomponent += 1

    key += '_{}'.format(id['gene'])
    try:
        if subj_cov > gene[key]:
            gene[key] = subj_cov
    except KeyError:
        gene[key] = subj_cov
        ngene += 1

    key += '_{}'.format(id['isoform'])
    try:
        if subj_cov > isoform[key]:
            isoform[key] = subj_cov
    except KeyError:
        isoform[key] = subj_cov
        nisoform += 1

    # if nhits >= 1000:
    #     break

print('{} blast hists processed'.format(nhits))
print('  {} clusters'.format(ncluster))
print('  {} components'.format(ncomponent))
print('  {} genes'.format(ngene))
print('  {} isoforms'.format(nisoform))

for key in cluster.keys():
    print('cluster:{}\t{:.3f}'.format(key, cluster[key]))

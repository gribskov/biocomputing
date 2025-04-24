"""=================================================================================================
Use blast results to identify taxonomic lineage for queries. This can be used to identify
contaminants in assemblies. Uses the JGI taxonomy servicelook at version before 4/24/2025 for
NCBI code (did not use because to excessive number of server errors)

blast_to_taxonomy.py <blast_file> <good_sequences>

Michael Gribskov     13 March 2025
================================================================================================="""
import sys
import json
from collections import defaultdict
from time import sleep
from urllib.parse import quote_plus

import requests
from lxml import etree


def send_query_jgi(tax, errors):
    """---------------------------------------------------------------------------------------------
    Deprecated. Not updated to query based on taxid
    send a single query to the jgi taxonomy server (usually you want to
    send in blocks, see
    send_query_block_jgi())

    :param tax: string      The stitle string from a blast search result
    :param errors: list     entries record queries that returned errors
    :return: string         full taxonomy string, semicolon delimited
    ---------------------------------------------------------------------------------------------"""
    url = 'https://taxonomy.jgi.doe.gov/name/simple/' + tax.strip().replace(' ', '_')

    response = None
    try:
        response = requests.post(url, {})
    except Exception as err:
        print(f'something goes wrong: {err}\n{response.text}')

    j = json.loads(response.text)
    taxstr = ''
    try:
        q = j[list(j)[0]]
        taxstr = ';'.join([q[l]['name'] for l in reversed(q.keys())
                           if l not in ('level', 'mononomial', 'name', 'tax_id')])

    except Exception as err:
        errors.append({'taxon': tax, 'error': err, 'response': response.text})

    return taxstr


def send_query_block_jgi(qlist, taxa):
    """---------------------------------------------------------------------------------------------
    send a list of sequences to the bjgi taxonomy services as a single post query

    :param qlist: list      taxonomic names from the stitle string from a blast search result
    :param taxa: dict       taxonomic information with blast taxa name as key; created by blast_read
    :return: int            number of queries successfully translated
    ---------------------------------------------------------------------------------------------"""
    url_taxstr = 'https://taxonomy.jgi.doe.gov/name/simple/'
    url_taxid = 'https://taxonomy.jgi.doe.gov/id/'
    block = ','.join(str(taxid) for taxid in qlist)
    query = url_taxid + block

    response = None
    j = None
    try:
        response = requests.post(query)
        j = json.loads(response.text)
    except Exception as err:
        print(f'something goes wrong: {err}\n{response.text}')

    ok_n = 0
    if len(j) != len(qlist):
        print(f'query and result are different lengths')
    oldtlen = len(taxa)
    i = 0
    for query in j:
        iquery =int(query)
        if qlist[i] != iquery:
            print(f'query/response mismatch')
        i += 1

        # try:
        taxinfo = j[query]
        if 'error' in taxinfo:
            # taxids that are not found in database
            taxa[iquery]['group'] = 'error'
            continue
        try:
            taxa[iquery]['lineage'] = ';'.join([taxinfo[l]['name'] for l in reversed(taxinfo.keys())
                           if l not in ('level', 'mononomial', 'name', 'tax_id')])
        except Exception as err:
            taxstr = 'Error;Taxon_not_found'

    return ok_n


def blast(blastfile):
    """---------------------------------------------------------------------------------------------
    generator that returns one blast result at each call
    format:
    # TRINITY_DN644_c0_g2_i1	3590	894	535	UniRef90_M1BIC5	388	269	388	120	100	624	4.21e-68	UniRef90_M1BIC5 Ninja-family protein n=2 Tax=Solanum TaxID=4107 RepID=M1BIC5_SOLTU
    #  qseqid qlen qstart qend sseqid slen sstart send length pident score evalue stitle
    :param blastfile:   string    path to blast search result
    :yield: dict        blast hit information for one hit
    ---------------------------------------------------------------------------------------------"""
    column = ['qid', 'qlen', 'qbegin', 'qend',
              'sid', 'slen', 'sbegin', 'send',
              'allen', 'pident', 'score', 'evalue', 'stitle']
    column_n = len(column)
    blast = open(blastfile, 'r')
    for line in blast:
        field = line.rstrip().split('\t', maxsplit=column_n - 1)
        yield {column[i]: field[i] for i in range(len(field))}

    blast.close()
    return


def blast_parse_stitle(stitle):
    """---------------------------------------------------------------------------------------------
    get the Tax string and numerical TaxID from the subject title field. This is specific for 
    information in the unir3ef databases
    
    UniRef90_A0A8T2XXL3 Vacuolar cation/proton exchanger n=2 Tax=Populus deltoides TaxID=3696 RepID=A0A8T2XXL3_POPDE
    ------------------- --------------------------------   -     -----------------       ----       ----------------
    urefid              description                        n     taxstr                  taxid      repid
    
    :param stitle: string       subject title from blast search
    :return: dict               keys: urefid, description, n, taxstr, taxid, repid
    ---------------------------------------------------------------------------------------------"""
    urefid, rest = stitle.split(' ', maxsplit=1)

    npos = rest.find('n=')
    taxpos = rest.find('Tax=')
    taxidpos = rest.find('TaxID')
    repidpos = rest.find('RepID')

    return {'urefid':      urefid,
            'description': rest[:npos - 1],
            'n':           int(rest[npos + 2:taxpos - 1]),
            'taxstr':      rest[taxpos + 4:taxidpos - 1],
            'taxid':       int(rest[taxidpos + 6:repidpos - 1]),
            'repid':       rest[repidpos + 6:]
            }


def blast_read(blastfile, trinity_level=3):
    """---------------------------------------------------------------------------------------------
    Read the taxonmic string from the blast search and make a list of unique taxa. Spaces in the
    taxon name are converted to underlines, after filtering parenthesized phrases

    :param blastfile: string    path to file with blast result
    :return: dict               taxid is the key, lineage is empty to be filled later
    ---------------------------------------------------------------------------------------------"""
    nline = 0
    blasthits = []
    taxa = defaultdict(lambda: {'taxstr': [], 'lineage': '', 'group': 'other', 'count': 0})

    for hit in blast(blastfile):
        nline += 1

        # save blast hit with truncated trinity ID and parsed information from stitle
        blasthits.append(hit)
        blasthits[-1]['truncated_id'] = trinity_filter_id(hit['qid'], trinity_level)

        info = blast_parse_stitle(hit['stitle'])
        for k in info:
            blasthits[-1][k] = info[k]

        taxid = info['taxid']
        taxstr = info['taxstr']

        taxa[taxid]['count'] += 1
        if taxstr not in taxa[taxid]['taxstr']:
            taxa[taxid]['taxstr'].append(taxstr)
        print(f'{nline:6d}\t{taxid}\t{taxstr}')

        # # debug
        # if nline > 10000:
        #     break

    print(f'\nblast results processed: {nline}')
    print(f'unique taxa {len(taxa)}')

    return taxa, blasthits


def blast_trinity_cluster(blasthits):
    """---------------------------------------------------------------------------------------------
    gather results for trinity ids clustered at a specified level. each truncated trinity id is a
    list of blast results (dicts) with the taxid extracted and added to the dicts

    :param blasthits: list          dict of each blasthit from the search, see blast_read()
    :return: dict                   key is truncated trinity ID, value is list of hits
    ---------------------------------------------------------------------------------------------"""
    cluster = defaultdict(lambda: {'group':'', 'member':[]})
    for hit in blasthits:
        cluster[hit['truncated_id']]['member'].append(hit)

    return cluster


def assign_groups(taxa, gorder, group):
    """---------------------------------------------------------------------------------------------
    use the search terms in group to assign assemblies to defined groups. sequences not matching
    are defined as other. error sequences are defined as 'error'

    :param taxa: dict       content: {'lineage': , 'group': , 'count': }
    :param group: dict      key is group name, value is taxonomic term to match
    :return: None
    ---------------------------------------------------------------------------------------------"""
    count = defaultdict(int)

    found = False
    g = None
    for t in taxa:
        this = taxa[t]
        for g in gorder:
            found = False
            if group[g] in this['lineage']:
                this['group'] = g
                found = True
                break
        if not found:
            this['group'] = 'other'

        count[g] += 1

    return count


def trinity_filter_id(tid, level=3, remove_trinity=True):
    """---------------------------------------------------------------------------------------------
    filter a trinity id to truncate it at the desired level, and optionally drop the redundant
    TRINITY_ prefix

    level 4: all parts, TRINITY_DN221212_c0_g1_i1
    level 3: gene,      TRINITY_DN221212_c0_g1
    level 2: component, TRINITY_DN221212_c0
    level 1: bundle,    TRINITY_DN221212

    :param tid: string                  trinity sequence ID
    :param level: int                   level to trim trinity IDs, default is gene level
    :param remove_trinity: boolean      remove the use TRINITY_ prefix if true
    :return: string                     trinity ID trimmed at the desired level
    ---------------------------------------------------------------------------------------------"""
    first = 0
    last = level + 1
    if remove_trinity:
        first = 1

    field = tid.split('_')
    return '_'.join(field[first:last])


def group_choose(obs_group):
    """---------------------------------------------------------------------------------------------
    assign the merged isoforms to a group based on the counts of the groups (obs_group).
    The rules are:
    1) assign as plant if there are any plants
    2) assign as highest count
    3) if it's a tie, assign as other

    :param merged: list of Cluster      annotated merged clusters
    :return: string                     count of number assigned to each group
    ---------------------------------------------------------------------------------------------"""
    gbest = ''
    n = 0
    
    for g in sorted(obs_group, key=lambda g: obs_group[g], reverse=True):
        if g == 'plant':
            return 'plant'

        if obs_group[g] == n:
            gbest = 'other'
        else:
            if obs_group[g] > n:
                gbest = g
                n = obs_group[g]

    return gbest


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    blastfile = sys.argv[1]
    goodfile = sys.argv[2]

    # define groups - each groups is defined by one taxonomic term, anything undefined is 'other'
    # it is important to process the most specific keywords first, gorder is the processing order
    group = {'error': 'Error', 'plant': 'Viridiplantae', 'virus': 'Viruses', 
            'arthropod': 'Arthropoda', 'fungi': 'Fungi', 'bacteria': 'Bacteria', 
            'opisthokonta':'Opisthokonta', 'sar':'Sar', 'discoba':'Discoba', 
            'amoebozoa': 'Amoebozoa', 'other':'Other'}
    gorder = ['error', 'plant', 'other', 'virus', 'arthropod', 'fungi', 'bacteria', 
            'opisthokonta', 'sar', 'discoba', 'amoebozoa', 'other']

    good_n = 0
    bad_n = 0
    unknown_n = 0
    nquery = 0

    taxa, blasthits = blast_read(blastfile)
    ntaxa = 0
    n_success = 0
    block = 250
    qlist = []
    for hit in taxa.keys():
        ntaxa += 1
        qlist.append(hit)

        if ntaxa % block:
            continue

        n_success += send_query_block_jgi(qlist, taxa)
        qlist = []

    if qlist:
        n_success += send_query_block_jgi(qlist, taxa)

    # assign to groups based on lineage
    groupcount = assign_groups(taxa, gorder, group)

    # ----------------------------------------------------------------------------------------------
    # results
    # ----------------------------------------------------------------------------------------------
    print(f'\nCounts per group found in {blastfile}')
    for g in groupcount:
        print(f'\t{g}\t{groupcount[g]}')

    good = open(goodfile, 'w')

    # species with top counts in blast result
    ntop = 100
    good.write(f'\n! top {ntop} taxa:\n')
    print(f'\ntop {ntop} taxa:')
    n = 0

    for t in sorted(taxa, key=lambda c: taxa[c]['count'], reverse=True):
        species = taxa[t]['lineage'].split(';')
        good.write(f"{n:3d}{taxa[t]['count']:10d}{taxa[t]['group']:>14s}  {species[-1]}\ttaxid={t:d}\n")
        print(f"{n:3d}{taxa[t]['count']:10d}{taxa[t]['group']:>14s}  {species[-1]}\ttaxid={t:d}")
        n += 1
        if n == ntop:
            break

    # assign groups to trinity clusters
    trinity_cluster = blast_trinity_cluster(blasthits)
    for tc in trinity_cluster:
        obs_group = defaultdict(int)
        for iso in trinity_cluster[tc]['member']:
            g = taxa[iso['taxid']]['group']
            obs_group[g] += 1
        trinity_cluster[tc]['group'] = group_choose(obs_group)

    # write out by groups
    good.write(f'\n! Taxa by group\n')
    group_old = ''
    nout = 0
    for c in sorted(trinity_cluster, key=lambda c: trinity_cluster[c]['group']):
        if trinity_cluster[c]['group'] != group_old:
            group_old = trinity_cluster[c]['group']
            good.write(f'\n! group:{group_old}\tsequences:{groupcount[group_old]}\n')
            print(f'\n! group:{group_old}\tsequences:{groupcount[group_old]}')

        good.write(f"!\t{c}\t{trinity_cluster[c]['group']}\t{len(trinity_cluster[c]['member'])}\n")
        print( f"!\t{c}\t{trinity_cluster[c]['group']}\t{len(trinity_cluster[c]['member'])}")
        for m in trinity_cluster[c]['member']:
            lineage = taxa[m['taxid']]['lineage']
            nout += 1
            print(f"{m['qid']}\t{trinity_cluster[c]['group']}\t{m['sid']}\t{m['evalue']}\t{m['taxstr']}\t{m['description']}")
            good.write(f"{m['qid']}\t{trinity_cluster[c]['group']}\t{m['sid']}\t{m['evalue']}\t{m['taxstr']}\t{m['description']} {lineage[:25]}\n")

    print(f'\n{nout} results written to {goodfile}')
    good.close()

    exit(0)

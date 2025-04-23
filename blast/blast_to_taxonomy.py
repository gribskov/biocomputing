"""=================================================================================================
Use blast results to identify taxonomic lineage for queries. This can be used to identify
contaminants in assemblies. Uses the JGI taxonomy service, but there is code here for NCBI as well.

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


def send_query(url, params, retry):
    """---------------------------------------------------------------------------------------------
    Send query and decode XML response
    :param url: string      URL of endpoint
    :param params: dict     parameters for query
    :param retry: list      list of failed queries
    :return: etree          parsed xml response
    ---------------------------------------------------------------------------------------------"""
    sleep(0.12)

    try:
        response = requests.post(url, params)
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        print("HTTP Error: {err}")
        print(err.args[0])
        retry.append(params)
        return None
    except requests.exceptions.ConnectionError as err:
        print(f"Connection error: {err}")
        retry.append(params)
        return None
    except requests.exceptions.Timeout as err:
        print(f"Timeout error: {err}")
        retry.append(params)
        return None
    except requests.exceptions.RequestException as err:
        print(f"request error: {err}")
        retry.append(params)
        return None

    if response.status_code != 200:
        print(f'something wrong')
        retry.append(params)
        return None

    try:
        print(f'Send_query status: {response.status_code}\t{response.text}')
    except Exception as err:
        print(f'send_query exception: {err}')
        return None

    # remove the encoding and doc type information before parsing with etree
    xmlstring = response.text[response.text.find('.dtd">') + 6:]
    xml = etree.fromstring(xmlstring)
    if xml is None:
        print(f'send_query: no xml document found {params}\n{xmlstring}')

    # check if search term not found
    notfound = xml.xpath('//PhraseNotFound')
    if len(notfound):
        print(f'not found: {xml}')
        retry.append(params)
        return None

    return xml


def send_query_jgi(tax, errors):
    """---------------------------------------------------------------------------------------------
    send a single query to the jgi taxonomy server (usually you want to send in blocks, see
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
    send a list of seqeunces to the bjgi taxonomy services as a single post query

    :param qlist: list        taxonomic names from the stitle string from a blast search result
    :param taxa: dict       taxonomic information with blast taxa name as key; created by blast_read
    :return: int            number of queries successfully translated
    ---------------------------------------------------------------------------------------------"""
    url = 'https://taxonomy.jgi.doe.gov/name/simple/'
    block = ','.join(qlist)
    query = url + block

    response = None
    j = None
    try:
        response = requests.post(query)
        j = json.loads(response.text)
    except Exception as err:
        print(f'something goes wrong: {err}\n{response.text}')

    ok_n = 0
    oldtlen = len(taxa)
    for query in j:

        try:
            taxinfo = j[query]
            taxstr = ';'.join([taxinfo[l]['name'] for l in reversed(taxinfo.keys())
                               if l not in ('level', 'mononomial', 'name', 'tax_id')])
        except TypeError:
            taxstr = 'Error;Taxon_not_found'

        try:
            taxstr = taxstr.replace(' ', '_')
            if not taxstr:
                taxstr = 'None'
            taxa[query]['lineage'] = taxstr
            if oldtlen != len(taxa):
                # some queries will trigger multiple responses, particularly ones with /
                # these queries are filtered in blast_read(), but there could be other
                # names that also do it
                print('Error:unknown_taxon attempting to extend taxa list. check filtering')
                oldtlen = len(taxa)
            ok_n += 1

        except Exception as err:
            # any other errors
            errstr = f'{err}'
            taxa[query]['lineage'] = f'Error; {errstr.replace(' ', '_')}'

    return ok_n


def get_taxinfo(xml):
    """---------------------------------------------------------------------------------------------
    Not used
    look up the target taxonomic name at the NCBI taxonomy database
    :param xml: etree       parsed xml response
    :return: int            NCBI TaxID
    ---------------------------------------------------------------------------------------------"""
    taxidlist = xml.xpath('//IdList/Id')

    taxid = None
    if taxidlist:
        taxid = int(taxidlist[0].text)
        # print(f'taxid:{taxid}')

    if not taxid:
        print(f'taxid:{taxid}\nxml:[{xml.text}]')
    return taxid


def get_webenv(xml):
    """---------------------------------------------------------------------------------------------
    Not used
    get the webenv ID and QueryKey from the xml

    :param xml: etree       parsed xml response
    :return: str, int       webenvID, QueryKey
    ---------------------------------------------------------------------------------------------"""
    webenv = xml.xpath('//WebEnv')[0].text
    querykey = int(xml.xpath('//QueryKey')[0].text)

    return webenv, querykey


def get_lineage(xml):
    """---------------------------------------------------------------------------------------------
    Not used
    get the lineage from the NCBI efetch xml result

    :param xml:
    :return:
    ---------------------------------------------------------------------------------------------"""
    lineage = 'unknown'

    try:
        lineage = xml.xpath('//Lineage')[0].text
    except Exception as err:
        print(f'{xml.text}\n{err}')

    return lineage


def get_lineage_from_taxonomydb(tax, retry):
    """---------------------------------------------------------------------------------------------
    Not used
    Give a taxonomy string such as Solanum commersonii or Solanoideae, get the taxonomic lineage
    from the NCBI taxonomy database

    This code works but the NCBI servers generate too many errors to use it for thousands of queries

    :param tax: string      taxonomic description
    :param retry: list      list for failed queries (search parameters are stored)
    :return: string         compleate lineage (comma separated)
    ---------------------------------------------------------------------------------------------"""
    api_key = 'f6646b2e3140e28ec200c579b03f21759107'
    ncbi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    efetch = 'efetch.fcgi'
    esearch = 'esearch.fcgi?'

    # query taxonomy db using taxonomic phrase from blast result
    params = {'api_key':    api_key,
              'db':         'taxonomy',
              'field':      'Organism',
              'term':       quote_plus(tax.strip()),
              'usehistory': 'y'}
    xml = send_query(ncbi + esearch, params, retry)

    # if xml is None, the query did not succeed, it is stored in retry, skip to next query
    if xml is None:
        return None

    if xml.findtext('PhraseNotFound'):
        print(f'1 Not found ({tax}). \n{xml}')

    taxid = get_taxinfo(xml)
    webenv, querykey = get_webenv(xml)

    # retrieve full information from taxonomy db

    params = {'api_key':    '645c588aed8ff727e6ed8059d10e7db2ea09',
              'db':         'taxonomy',
              'TaxID':      taxid,
              'usehistory': 'y',
              'WebEnv':     webenv,
              'query_key':  querykey}
    xml = send_query(ncbi + efetch, params, retry)
    if xml is None:
        return None

    lineage = get_lineage(xml)
    if lineage == 'unknown':
        print(f'taxid:{taxid}\nxml:|{xml.text}|')

    return lineage


def blast(blastfile):
    """---------------------------------------------------------------------------------------------
    generator that returns one blast result at each call
    format

    # TRINITY_DN644_c0_g2_i1	3590	894	535	UniRef90_M1BIC5	388	269	388	120	100	624	4.21e-68	UniRef90_M1BIC5 Ninja-family protein n=2 Tax=Solanum TaxID=4107 RepID=M1BIC5_SOLTU
    #  qseqid qlen qstart qend sseqid slen sstart send length pident score evalue stitle
    :param blastfile: string    path to blast search result
    :return:
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


def blast_get_tax(stitle):
    """---------------------------------------------------------------------------------------------
    get the Tax string and numerical ID from the subject title field

    :param stitle: string       subject title from blast search
    :return: string, int        Tax string, taxid
    ---------------------------------------------------------------------------------------------"""
    # taxid is in a qualifier 'TaxID=xxx '
    taxidpos = stitle.find('TaxID') - 1
    taxid_end = stitle.find('RepID') - 1
    taxid = int(stitle[taxidpos + 7:taxid_end])

    taxpos = stitle.find('Tax=') + 4
    return stitle[taxpos:taxidpos], taxid


def blast_read(blastfile):
    """---------------------------------------------------------------------------------------------
    Read the taxonmic string from the blast search and make a list of unique taxa. Spaces in the
    taxon name are converted to underlines, after filtering parenthesized phrases

    :param blastfile: string    path to file with blast result
    :return: dict               taxonomic string is the key, lineage is empty to be filled later
    ---------------------------------------------------------------------------------------------"""
    nline = 0
    taxa = defaultdict(lambda: {'taxid': 0, 'lineage': '', 'group': 'other', 'count': 0})

    for hit in blast(blastfile):
        nline += 1
        tax, taxid = blast_get_tax(hit['stitle'])

        # trim some unecessary/confusing parts from the taxonomy string
        if tax.find('(') != -1:
            # trim off information in parentheses
            # this information is usually not part of the actual taxonomic name
            tax = tax[:tax.find('(') - 1]
            # print(f'\t=> trimmed: {tax}')

        # convert spaces to underline
        tax = tax.strip().replace(' ', '_')
        tax = tax.replace('/', '_')

        # remove strain information (sometimes causes multiple returns from taxonomy server)
        strain = tax.find('_str.')
        if strain > -1:
            tax = tax[:strain]
        # sp = tax.find('_sp.')
        # if sp > -1:
        #     tax = tax[:sp]
        uncultured = tax.find('uncultured_')
        if uncultured > -1:
            tax = tax[uncultured + 11:]
        unclassifed = tax.find('unclassified_')
        if unclassifed > -1:
            tax = tax[unclassifed + 13:]
        # if tax.find('str.')>-1:
        #     print(tax)
        taxa[tax]['count'] += 1
        taxa[tax]['taxid'] = taxid
        print(f'{nline:6d}\t{tax}')

        # debug
        if nline > 10000:
            break

    print(f'\nblast results processed: {nline}')
    print(f'unique taxa {len(taxa)}')

    return taxa


def blast_trinity_cluster(blastfile):
    """---------------------------------------------------------------------------------------------
    gather results for trinity ids clustered at a specified level. each truncated trinity id is a
    list of blast results (dicts) with the taxid extracted and added to the dicts

    :param blastfile: string        path to blast result file
    :return: dict                   key is truncated trinity ID, value is list of hits
    ---------------------------------------------------------------------------------------------"""
    cluster = defaultdict(lambda: [])
    for hit in blast(blastfile):
        tid = trinity_filter_id(hit['qid'], 3, True)
        taxstr, taxid = blast_get_tax(hit['stitle'])
        hit['taxid'] = taxid
        cluster[tid].append(hit)

    return cluster


def assign_groups(taxa, group):
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
        for g in group:
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
    for iso in merged:
        if iso['group'] == 'plant':
            count['plant'] += 1
            continue

        # gsort = sorted(this.ngroup, key=lambda g: this.ngroup[g], reverse=True)
        # if len(gsort) == 1:
        #     gbest = gsort[0]
        # elif gsort[0] > gsort[1]:
        #     gbest = gsort[0]
        # else:
        #     gbest = 'other'
        #
        # this.group = gbest
        # count[gbest] += 1

    return count

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
            return 'other'
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
    group = {'error':     'Error', 'plant': 'Viridiplantae', 'virus': 'Viruses',
             'arthropod': 'Arthropoda', 'fungi': 'Fungi', 'bacteria': 'Bacteria'}
    gorder = {'error': 1, 'plant': 3, 'virus': 4, 'arthropod': 5, 'fungi': 6, 'bacteria': 7, 'other': 2}

    good_n = 0
    bad_n = 0
    unknown_n = 0
    nquery = 0

    taxa = blast_read(blastfile)
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
    groupcount = assign_groups(taxa, group)
    taxa_invert = {}
    for t in taxa:
        taxa_invert[taxa[t]['taxid']] = taxa[t]
        taxa_invert[taxa[t]['taxid']]['taxstr'] = t
        if t.find('Glarea')>-1 or taxa[t]['taxid']==1104152 or taxa[t]['taxid']==101852:
            print(f'{taxa[t]}')

    for t in sorted(taxa_invert):
        if t == 1104152:
            print(f'{t}\t{taxa_invert[t]}')

    # ----------------------------------------------------------------------------------------------
    # results
    # ----------------------------------------------------------------------------------------------
    # print(f'\nCounts per group found in {blastfile}')
    # for g in groupcount:
    #     print(f'\t{g}\t{groupcount[g]}')
    #
    # good = open(goodfile, 'w')
    #
    # # species with top counts in blast result
    # ntop = 100
    # good.write(f'\n! top {ntop} taxa:\n')
    # print(f'\ntop {ntop} taxa:')
    # n = 0
    #
    # for t in sorted(taxa, key=lambda c: taxa[c]['count'], reverse=True):
    #     good.write(f'!{n:3d}{t:>40s}{taxa[t]['group']:>10s}{taxa[t]['count']:>10d}  {taxa[t]['lineage']}\n')
    #     print(f'{n:3d}{t:>40s}{taxa[t]['group']:>10s}{taxa[t]['count']:10d}  {taxa[t]['lineage']}')
    #     n += 1
    #     if n == ntop:
    #         break
    #
    # good.write(f'\n! Taxa by group\n')
    # group_old = ''
    # n = 0
    # for tax in sorted(taxa, key=lambda t: (gorder[taxa[t]['group']],taxa[t]['lineage'])):
    #     n += 1
    #     this = taxa[tax]
    #     good.write(f"{n:5d}{this['group']:>12s} {tax:>68s}\t{this['lineage']}\n")
    #
    #
    # print(f'\nresults written to {goodfile}')
    # good.close()

    # reprocess the blast results, pooling at trinity level and summarizing taxonomic information
    trinity_cluster = blast_trinity_cluster(blastfile)
    for tc in trinity_cluster:
        obs_group = defaultdict(int)
        for iso in trinity_cluster[tc]:
            g = taxa_invert[iso['taxid']]['group']
            obs_group[g] += 1
        best_group = group_choose(obs_group)




    exit(0)

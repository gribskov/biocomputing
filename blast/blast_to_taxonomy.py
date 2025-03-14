"""=================================================================================================
Use blast results to identify taxonomic lineage for queries. This can be used to identify
contaminants in assemblies

<TaxaSet>
<Taxon>
<TaxId>12183</TaxId>
<ScientificName>Potato virus X</ScientificName>
<OtherNames>
<Acronym>PVX</Acronym>
<EquivalentName>potato virus X, PVX</EquivalentName>
</OtherNames>
<ParentTaxId>12176</ParentTaxId>
<Rank>species</Rank>
<Division>Viruses</Division>
<GeneticCode>
<GCId>1</GCId>
<GCName>Standard</GCName>
</GeneticCode>
<MitoGeneticCode>
<MGCId>0</MGCId>
<MGCName>Unspecified</MGCName>
</MitoGeneticCode>
<Lineage>Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Alsuviricetes; Tymovirales; Alphaflexiviridae; Potexvirus</Lineage>
<LineageEx>
<Taxon>
<TaxId>10239</TaxId>
<ScientificName>Viruses</ScientificName>
<Rank>no_rank</Rank>
</Taxon>
<Taxon>
<TaxId>2559587</TaxId>
<ScientificName>Riboviria</ScientificName>
<Rank>realm</Rank>
</Taxon>
<Taxon>
<TaxId>2732396</TaxId>
<ScientificName>Orthornavirae</ScientificName>
<Rank>kingdom</Rank>
</Taxon>
<Taxon>
<TaxId>2732406</TaxId>
<ScientificName>Kitrinoviricota</ScientificName>
<Rank>phylum</Rank>
</Taxon>
<Taxon>
<TaxId>2732461</TaxId>
<ScientificName>Alsuviricetes</ScientificName>
<Rank>class</Rank>
</Taxon>
<Taxon>
<TaxId>675063</TaxId>
<ScientificName>Tymovirales</ScientificName>
<Rank>order</Rank>
</Taxon>
<Taxon>
<TaxId>675064</TaxId>
<ScientificName>Alphaflexiviridae</ScientificName>
<Rank>family</Rank>
</Taxon>
<Taxon>
<TaxId>12176</TaxId>
<ScientificName>Potexvirus</ScientificName>
<Rank>genus</Rank>
</Taxon>
</LineageEx>
<CreateDate>1995/02/27 09:24:00</CreateDate>
<UpdateDate>2020/04/07 15:24:16</UpdateDate>
<PubDate>1993/04/23 01:00:00</PubDate>
</Taxon>
</TaxaSet>

Michael Gribskov     13 March 2025
================================================================================================="""

import requests
from time import sleep
from urllib.parse import quote_plus
from lxml import etree
from blast import Blast


def send_query(url, params, retry):
    """---------------------------------------------------------------------------------------------
    Send query and decode XML response
    :param url: string      URL of endpoint
    :param params: dict     parameters for query
    :return: etree          parsed xml response
    ---------------------------------------------------------------------------------------------"""
    sleep(0.12)
    # try:
    #     r = requests.get(url, timeout=1)
    #     r.raise_for_status()
    # except requests.exceptions.HTTPError as errh:
    #     print("HTTP Error")
    #     print(errh.args[0])

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


def get_taxinfo(xml):
    """---------------------------------------------------------------------------------------------
    look up the target taxonomic name at the NCBI taxonomy database
    :param xml: etree       parsed xml response
    :return: int            NCBI TaxID
    ---------------------------------------------------------------------------------------------"""
    taxidlist = xml.xpath('//IdList/Id')

    if taxidlist:
        taxid = int(taxidlist[0].text)
        # print(f'taxid:{taxid}')

    if not taxid:
        print(f'taxid:{taxid}\nxml:[{xml.text}]')
    return taxid


def get_webenv(xml):
    """---------------------------------------------------------------------------------------------
    get the webenv ID and QueryKey from the xml

    :param xml: etree       parsed xml response
    :return: str, int       webenvID, QueryKey
    ---------------------------------------------------------------------------------------------"""
    webenv = querykey = None
    webenv = xml.xpath('//WebEnv')[0].text
    querykey = int(xml.xpath('//QueryKey')[0].text)

    return webenv, querykey


def get_lineage(xml):
    """---------------------------------------------------------------------------------------------
    get the lineage from the efetch xml result

    :param xml:
    :return:
    ---------------------------------------------------------------------------------------------"""
    lineage = 'unknown'

    try:
        lineage = xml.xpath('//Lineage')[0].text
    except Exception as err:
        print(f'{xml.text}\n{err}')

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
    get the Tax string from the subject title field

    :param stitle: string       subject title from blast search
    :return: string             Tax string
    ---------------------------------------------------------------------------------------------"""
    taxpos = stitle.find('Tax=') + 4
    taxidpos = stitle.find('TaxID')
    return stitle[taxpos:taxidpos]


def get_lineage_from_taxonomydb(tax, retry):
    """---------------------------------------------------------------------------------------------
    Give a taxonomy string such as Solanum commersonii or Solanoideae, get the taxonomic lineage
    from the NCBI taxonomy database

    :param tax: string      taxonomic description
    :param retry: list      list for failed queries (search parameters are stored)
    :return: string         compleate lineage (comma separated)
    ---------------------------------------------------------------------------------------------"""
    api_key = 'f6646b2e3140e28ec200c579b03f21759107'
    ncbi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    efetch = 'efetch.fcgi'
    esearch = 'esearch.fcgi?'

    # query taxonomy db using taxonomic phrase from blast result
    params = {'api_key':    '645c588aed8ff727e6ed8059d10e7db2ea09',
              'db':         'taxonomy',
              'term':       quote_plus(tax),
              'usehistory': 'y'}
    xml = send_query(ncbi + esearch, params, retry)

    # if xml is None, the query did not succeed, it is stored in retry, skip to next query
    # TODO check if timeout can be extended
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


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    blastfile = 'data/c16c31.trinity_uniref_1e-10.dmndblastx'
    taxa = {}
    good_n = 0
    bad_n = 0
    retry = []
    nline = 0
    nquery = 0
    for hit in blast(blastfile):
        nline += 1
        tax = blast_get_tax(hit['stitle'])
        print(f'{nline}\t{tax}')
        if tax.find('(') != -1:
            tax = tax[:tax.find('(') - 1]
            print(f'\t=> trimmed: {tax}')

        if tax in taxa:
            # skip if lineage is already known
            continue

        nquery += 1
        lineage = get_lineage_from_taxonomydb(tax, retry)
        if lineage:
            print(f'\t{nquery}\t{lineage[:60]}')
        else:
            print(f'\t{nquery}\tNone')
        if lineage is None:
            # query failed, go on to next line
            continue

        # print(f'tax:{tax}| {lineage}')
        if lineage is None:
            print(f'no such taxon: {tax}')
            print(f'{hit}')
            continue

        elif lineage == 'unknown':
            print(f'unknown taxon: {tax}\n{hit}')
            continue

        # plants are good, everything else is bad
        try:
            if lineage.index('Viridiplantae'):
                taxa[tax] = {'good': 1, 'lineage': lineage}
                good_n += 1
        except ValueError:
            taxa[tax] = {'good': 0, 'lineage': lineage}
            bad_n += 1

    print(f'good taxa: {good_n}\bbad taxa: {bad_n} found in {blastfile}')

    for tax in sorted(taxa, key=lambda t: taxa[t]['lineage']):
        if tax['good']:
            print(f'tax:{tax}')

    exit(0)

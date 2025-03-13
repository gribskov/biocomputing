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
from urllib.parse import quote_plus
from lxml import etree


def send_query(url, params):
    """---------------------------------------------------------------------------------------------
    Send query and decode XML response
    :param url: string      URL of endpoint
    :param params: dict     parameters for query
    :return: etree          parsed xml response
    ---------------------------------------------------------------------------------------------"""
    response = requests.post(url, params)
    # print(f'\n post response text: {response.text}')
    # remove the encoding and doc type information before parsing with etree
    return etree.fromstring(response.text[response.text.find('.dtd">') + 7:])


def get_taxid(xml):
    """---------------------------------------------------------------------------------------------
    look up the target taxonomic name at the NCBI taxonomy database
    :param xml: etree       parsed xml response
    :return: int            NCBI TaxID
    ---------------------------------------------------------------------------------------------"""
    taxid = None
    taxidlist = xml.xpath('//IdList/Id')

    if taxidlist:
        taxid = int(taxidlist[0].text)
        # print(f'taxid:{taxid}')

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
    return xml.xpath('//Lineage')[0].text


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    api_key = '645c588aed8ff727e6ed8059d10e7db2ea09'
    ncbi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    efetch = 'efetch.fcgi'
    esearch = 'esearch.fcgi?'

    # query taxonomy db using taxonomic phrase from blast result
    target = "solanum"
    params = {'api_key':    '645c588aed8ff727e6ed8059d10e7db2ea09',
              'db':         'taxonomy',
              'term':       quote_plus(target),
              'usehistory': 'y'}
    xml = send_query(ncbi + esearch, params)

    taxid = get_taxid(xml)
    webenv, querykey = get_webenv(xml)

    # retrieve full information from taxonomy db

    params = {'api_key':    '645c588aed8ff727e6ed8059d10e7db2ea09',
              'db':         'taxonomy',
              'TaxID':      taxid,
              'usehistory': 'y',
              'WebEnv':     webenv,
              'query_key':  querykey}
    xml = send_query(ncbi + efetch, params)
    lineage = get_lineage(xml)
    print(f'{target} => {lineage}')

    exit(0)

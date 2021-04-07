from bs4 import BeautifulSoup, Comment

# check response to original submission
doc = open('result.test.xml', 'r').read()
# doc = doc.replace('\\n', '')
soup = BeautifulSoup(doc, "xml")
# print(soup.prettify())

hits = soup.find_all('Hit')

blasthits = []
for hit in hits:

    this_hit = {'id':       hit.Hit_id.text, 'def':hit.Hit_def.text,
                'accession':hit.Hit_accession.text, 'len':hit.Hit_len.text,
                'hsp':      []}
    blasthits.append(this_hit)
    for hsp in hit.find_all('Hsp'):
        hsp_num = hsp.find('Hsp_num').text
        print('\n{} hsp num:{}'.format(hit.Hit_id.text, hsp_num))
        this_hsp = {'bitscore': hsp.find('Hsp_bit-score').text,
                    ':hspscore':hsp.find('Hsp_score').text,
                    'evalue':   hsp.find('Hsp_evalue').text,
                    'identity': hsp.find('Hsp_identity').text,
                    'positive': hsp.find('Hsp_positive').text,
                    'gaps':     hsp.find('Hsp_gaps').text,
                    'align_len':hsp.find('Hsp_align-len').text,
                    'qbegin':   hsp.find('Hsp_query-from').text,
                    'qend':     hsp.find('Hsp_query-to').text,
                    'sbegin':   hsp.find('Hsp_hit-from').text,
                    'send':     hsp.find('Hsp_hit-to').text}
        this_hit['hsp'].append(this_hsp)

exit(0)

"""=================================================================================================
pull the gene ontology GO terms for each uniref sequence out of the xml file
================================================================================================="""
import sys
import yaml
import json
from lxml import etree

# from io import StringIO

xfile = sys.argv[1]
# print('xml file {}'.format(xfile))
# uniref = etree.parse(xfile)
uniref = None
try:
    uniref = open(xfile, 'r')
except OSError:
    sys.stderr.write('Error opening xml file ({})'.format(xfile))
    exit(1)

# entries = blast.xpath('//entry')

nseq = 0
seq = {}
for line in uniref:
    line = line.strip()

    if line.startswith('<entry '):
        nseq += 1
        if nseq > 100:
            break
        field = line.split(' ')
        id = field[1].replace('id="', '').replace('"', '')
        print('{} {}'.format(nseq, id))
        seq[id] = []

    elif (line.startswith('<property type="GO Molecular Function')):
        value = line[line.index("value") + 10:]
        quote = value.index('"')
        field = line.split(' ')
        go = field[2].replace('value="', '').replace('"', '')
        print('\tMF {}'.format(value[:quote]))
        seq[id].append('MF:{}'.format(value[:quote]))
    elif (line.startswith('<property type="GO Biological Process')):
        value = line[line.index("value") + 10:]
        quote = value.index('"')
        field = line.split(' ')
        go = field[2].replace('value="', '').replace('"', '')
        print('\tBP {}'.format(value[:quote]))
        seq[id].append('BP:{}'.format(value[:quote]))
    elif (line.startswith('<property type="GO Cellular Component')):
        value = line[line.index("value") + 10:]
        quote = value.index('"')
        field = line.split(' ')
        go = field[2].replace('value="', '').replace('"', '')
        print('\tCC {}'.format(value[:quote]))
        seq[id].append('CC:{}'.format(value[:quote]))

for id in seq:
    print(id, seq[id])

sys.stdout.write(json.dump(seq))
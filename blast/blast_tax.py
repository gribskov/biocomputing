"""=================================================================================================
analyze taxonomy of blast hits to identify contaminants

Michael Gribskov     18 February 2025
================================================================================================="""
import sys
from blast import Blast
import re
from collections import defaultdict


def get_tax_and_id(string):
    taxpos = string.find('Tax=')
    taxidpos = string.find('TaxID=')
    repidpos = string.find('RepID=')

    tax = string[taxpos+4:taxidpos-1]
    taxid = string[taxidpos+6:repidpos-1]

    return tax, taxid


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    infile = sys.argv[1]
    sys.stdout.write('Blast search: {}\n\n'.format(infile))
    blast = Blast(file=sys.argv[1])

    taxre = re.compile(r'Tax=(.*) TaxID=(.*) ')
    sidre = re.compile(r'UniRef90_')
    stopword = ['UniRef90_[^ ]+', r'n=\d+', r'sp\.* (\d+)*', 'protein', 'uncharacterized',
                'family', '-*domain-containing', r'\(?fragment\)?', r'\(strain[^)]*\)',
                'Tax=', 'TaxID=', 'RepID=[^ ]+']
    stopre = re.compile('|'.join(stopword), re.I)
    # regex for splitID()
    idre = re.compile(r'>*TRINITY_DN([^_]+)_c(\d+)_g(\d+)_i(\d+)')

    fmt = 'qid qlen qbegin qend sid slen sbegin send alignlen pid score evalue stitle'
    # fmt = 'qname sname id alignlen mismatch gapopen qbeg qend sbeg send evalue bit_score'
    nfields = blast.setFormat(fmt)

    n = 0
    taxa = defaultdict(lambda: 1)
    while blast.next():
        n += 1
        tax, taxid = get_tax_and_id(blast.stitle)
        taxa[f'{tax}|{taxid}'] += 1
        # print(f'{blast.stitle}|{tax}|{taxid}|')
        # if n == 50000:
        #     break

    for t in sorted(taxa, key=lambda t: taxa[t]):
        print(f'{taxa[t]}\t{t}')

    exit(0)

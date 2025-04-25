"""=================================================================================================
check for duplicate names in a fasta file. if sequences are the same length keep the first one,
otherwise rename the second one and keep both.

Michael Gribskov     25 April 2025
================================================================================================="""
import sys

def getfasta(fastainfile):
    """---------------------------------------------------------------------------------------------
    generator to return successive fasta sequences
    :param fastainfile:
    :return:
    ---------------------------------------------------------------------------------------------"""
    fastain = open(fastainfile, 'r')

    oldid = ''
    entry = None
    for line in fastain:
        if line.startswith('>'):
            if oldid:
                yield entry
                oldid = ''
            entry = {}
            try:
                id, doc = line.rstrip().split(maxsplit=1)
            except ValueError:
                id = line.rstrip()
                doc = ''
            entry['id'] = id.replace('>', '')
            entry['len'] = 0
            entry['seq'] = ''
            oldid = entry['id']
        else:
            entry['seq'] += line.rstrip()
            entry['len'] = len(entry['seq'])

    if oldid:
        yield entry

    fastain.close()
    return None

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fastainfile = sys.argv[1]
    fastaoutfile = sys.argv[2]
    fastaout = open(fastaoutfile, 'w')

    for f in getfasta(fastainfile):
        print(f)


    fastaout.close()
    exit(0)

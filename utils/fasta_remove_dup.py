"""=================================================================================================
check for duplicate names in a fasta file. if sequences are the same length keep the first one,
otherwise rename the second one and keep both.

Michael Gribskov     25 April 2025
================================================================================================="""
import sys

from sympy.solvers.diophantine.diophantine import length


def getfasta(fastainfile):
    """---------------------------------------------------------------------------------------------
    generator to return successive fasta sequences
    :param fastainfile:
    :yield: dict            keys: id, doc, seq, len
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
            entry['doc'] = doc
            entry['len'] = 0
            entry['seq'] = ''
            oldid = entry['id']
        else:
            entry['seq'] += line
            entry['len'] = len(entry['seq']) - 1

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

    seqs = {}
    nread = 0
    nskipped = 0
    nwritten = 0
    for f in getfasta(fastainfile):
        nread += 1

        if f['id'] in seqs:
            # sequence has been seen before
            if seqs[f['id']] == f['len']:
                # skip sequences that are the same length
                nskipped += 1
                continue

            # create new non overlapping name
            n = 1
            id = f['id'] + f'.{n}'
            while id in seqs:
                n += 1
                id = f['id'] + f'.{n}'
            fastaout.write(f">{id} {f['doc']}\n{f['seq']}")
            seqs['id'] = f['len']

        else:
            fastaout.write(f">{f['id']} {f['doc']}\n{f['seq']}")
            seqs[f['id']] = f['len']

        nwritten += 1

    print(f'sequences read from {fastainfile}: {nread}')
    print(f'duplicate sequences: {nskipped}')
    print(f'sequences written to {fastaoutfile}: {nwritten}')

    fastaout.close()
    exit(0)

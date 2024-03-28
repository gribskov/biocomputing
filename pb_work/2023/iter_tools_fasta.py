from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a multiple sequence fasta file. yield header and sequence of each
    """
    fh = open(fasta_name)

    # the groupby function yields pairs of values, a boolean indicating whether
    # the line is a docline, and a list like object with the lines that are
    # either doclines (one line) or sequence (one or more lines)
    linetype = groupby(fh, lambda line: line[0] == ">")

    doc = ''
    seq = ''
    for isdocline, text in linetype:

        if isdocline:
            doc = list(text)[0][1:-1]
        else:
            seq += ''.join(list(text))
            seq.replace('\n','')
            yield doc, seq
            doc = ''
            seq = ''

fasta_name = '../../data/Trinity.fasta'
iter = fasta_iter(fasta_name)

for doc, seq in iter:
    print(f'{doc}\n{seq}')
exit(0)
"""=================================================================================================
select sequences from a fasta file based on a table of desired IDs

Michael Gribskov     10 March 2025
================================================================================================="""


def getfasta(fastafile):
    """---------------------------------------------------------------------------------------------
    Generator that returns the next fasta sequence from the fastafile. Formatiing of original file
    is retained - newlines and > are kept from original.

    :param fastafile: string    path to fasta sequence file
    :return: dict               keys: id, documentation, sequence
    ---------------------------------------------------------------------------------------------"""
    fasta = open(fastafile, 'r')
    entry = {'id': '', 'documentation': '', 'sequence': ''}
    for line in fasta:
        if line.startswith('>'):
            # id/documentation line
            if entry['sequence']:
                yield entry
                entry = entry = {'id': '', 'documentation': '', 'sequence': ''}
            field = line.split(maxsplit=1)  # newline not removed
            if len(field) > 1:
                entry['documentation'] = field[1]
            entry['id'] = field[0].replace('>','')

        else:
            # sequence line (retain newlines)
            entry['sequence'] += line

    if entry['sequence']:
        yield entry

    fasta.close()
    return


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    idfile = 'data/count_plus_tax.selected.list'
    column = 0
    fastafile = 'data/C16C31.trinity.fasta'
    outfile = 'data/C16C31.taxselected.fasta'

    # read in list of IDs
    id = open(idfile, 'r')
    ids = []
    id_n = 0
    for line in id:
        field = line.rstrip().split()
        ids.append(field[0])
        id_n += 1

    print(f'ids: {id_n} IDs read from {idfile}')
    id.close()

    # read fasta and write to output
    out = open(outfile, 'w')
    fasta_n = 0
    fasta_sel = 0
    for fasta in getfasta(fastafile):
        fasta_n += 1
        if not fasta_n % 1000:
            print('.', end='')
        if not fasta_n % 40000:
            print(f' {fasta_n}\t{fasta_sel}')
        if fasta['id'] in ids:
            fasta_sel += 1
            out.write(f">{fasta['id']} {fasta['documentation']}")
            out.write(f"{fasta['sequence']}")

    print(f'\nfasta sequences read:{fasta_n} from {fastafile}')
    print(f'fasta sequences written:{fasta_sel} to {outfile}')
    out.close()

    exit(0)

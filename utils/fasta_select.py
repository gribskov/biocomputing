"""-------------------------------------------------------------------------------------------------
select sequences from a list of IDs

Input id.list is a list of names, only the first token is used
Output Fasta format, 00 letters/line

usage
    fasta_select id.list *.fasta > selected.fasta

2 September 2018    Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import glob
import sys
import re
from sequence.fasta import Fasta


def read_id(idfile):
    """---------------------------------------------------------------------------------------------
    Read a file of IDs and store in a list.  Only the first token is used
    Lines beginning ! or # are skipped

    :param idfile: open filehandle
    :return: list of strings
    ---------------------------------------------------------------------------------------------"""
    nid = 0
    idlist = []
    for line in idfile:
        if line.startswith(('#', '!')):
            continue
        if line.rstrip():
            id = line.replace('"','').split()[0]
            if id == 'baseMean':
                continue
            if id not in idlist:
                idlist.append(id)

    return idlist


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    sys.stderr.write('fasta.select\n')

    linelen = 100
    trim = re.compile('path=.*')

    # idlist
    idlistname = sys.argv[1]
    sys.stderr.write('\tID list: {}\n'.format(idlistname))
    try:
        idlist = open(idlistname, 'r')
    except:
        sys.stderr.write('\nUnable to open ID list ({})\n'.format(idlistname))
        exit(1)

    #  default target file name
    target = '*.fasta'
    if len(sys.argv) > 2:
        target = sys.argv[2]
    sys.stderr.write('\ttarget file: {}\n\n'.format(target))

    # read the ID list
    id = read_id(idlist)
    sys.stderr.write('{} sequence IDs read from {}\n'.format(len(id), idlistname))

    # read the sequences and store all that match the IDs
    # duplicates will be stored twice
    n_match = {}        # per file number of sequences in list
    n_notmatch = {}     # per file number of sequences not in list
    n_sequence = {}     # per file number of sequences
    n_found = {}        # per ID, number of times found in all files
    n_file = 0
    n_total = 0
    for fastafile in glob.glob(target):
        fasta = Fasta()
        fasta.open(fastafile)

        n_sequence[fastafile] = 0
        n_match[fastafile] = 0
        n_notmatch[fastafile] = 0
        n_file += 1
        n_written = 0

        while fasta.next():
            n_sequence[fastafile] += 1
            n_total += 1

            if fasta.id in id:
                # desired selected sequences
                fasta.trimDocByRegex(trim)
                sys.stdout.write('{}\n'.format(fasta.format(linelen=100)))
                n_written += 1
                n_match[fastafile] += 1
                if fasta.id in n_found:
                    n_found[fasta.id] += 1
                else:
                    n_found[fasta.id] = 1

            else:
                # not selected sequence
                n_notmatch[fastafile] += 1

    sys.stderr.write('files read: {}\n'.format(n_file))
    sys.stderr.write('total sequences read: {}\n'.format(n_total))
    sys.stderr.write('total sequences written: {}\n'.format(n_written))

    sys.stderr.write('\nPer file\n')
    for fastafile in n_sequence:
        sys.stderr.write('{}\n'.format(fastafile))
        sys.stderr.write('\tsequences: {}\n'.format(n_sequence[fastafile]))
        sys.stderr.write('\tsequences matched: {}\n'.format(n_match[fastafile]))

exit(0)

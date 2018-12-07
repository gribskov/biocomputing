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
            id = line.split()[0]
            if id not in idlist:
                idlist.append(id)

    exit(0)
    return idlist


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    sys.stderr.write('fasta.select\n')

    linelen = 100

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
    n_match = {}
    n_unique = {}
    n_sequence = {}
    n_file = 0
    n_total = 0
    for fastafile in glob.glob(target):
        fasta = Fasta()
        fasta.open(fastafile)
        n_sequence[fastafile] = 0
        n_unique[fastafile] = 0
        n_match[fastafile] = 0
        n_file += 1

        while fasta.next():
            n_sequence[fastafile] += 1
            if fasta.id in id:
                n_match[fastafile] += 1
                sys.stdout.write('{}\n'.format(fasta.format(linelen=100)))
            else:
                n_uniqueperfile += 1


        n_total += n_perfile
        n_unique_total += n_uniqueperfile
        sys.stderr.write(
            '{}\t{}\t{}\t{}\t{}\t{}\n'.format(n_file, fastafile, n_perfile, n_uniqueperfile,
                                              n_total,
                                              n_unique_total))

    # write out sequences
    for seq in unique_seq:
        sys.stdout.write(unique_seq[seq])

exit(0)

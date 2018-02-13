"""=================================================================================================
Fasta sequence class.  Supports iteration over a multi-fasta file.  The object has the following
externally accessible attributes:

    filename    input filename
    fh          input filehandle
    id          ID string (without>)
    doc         documentation string
    seq         sequence (one letter code, sequence characters only)

Synopsis
    from fasta import Fasta

    fasta = Fasta()
    fasta.open('filename')
    while True:
        fasta.next():
        print(fasta.format(linelen=60))

TODO: Don't forget to put your name and the date here
================================================================================================="""


class Fasta:

# TODO: write the Fasta class here. You must provide
#   attributes: id, doc, seq
#   methods:
#       open(<filename string>)
#       next()                      reads the next seq internally
#       format(<linelen>)           format the sequence with <linelen> sequence characters/line
#                                   <linelen> is optional, default should be 60

# note that the sequences should be cleaned of spaces, Ns and *s

# =================================================================================================
# main program - do not change the code below
# ==================================================================================================
if __name__ == '__main__':

    outfile = 'hw2.out'
    try:
        out = open(outfile, 'w')
    except:
        print('unable to open {}'.format(outfile))
        exit(1)

    fasta = Fasta()
    fasta.open('short_broke.fa')

    n = 0
    while True:
        n += 1
        if fasta.next():
            print('sequence {}:{} length={}'.format(n, fasta.id, len(fasta.seq)))
            out.write('{}'.format(fasta.format(linelen=50)))
        else:
            break

        if n >= 10:
            break

    out.close()
    exit(0)

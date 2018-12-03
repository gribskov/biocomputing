"""-------------------------------------------------------------------------------------------------
deg_add_blast

Add blast search information to a tab delimited gene expression table.  The first column of the the
deg table must match the query name of the blast search

3 December 2018     Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys
from blast.blast import Blast

# ==================================================================================================
if __name__ == '__main__':

    blastfile = sys.argv[1]
    sys.stderr.write('deg_add_blast')
    sys.stderr.write('\tBlast file: {}\n'.format(blastfile))

    blast = Blast()
    blast.new(blastfile)

    # format (diamond)

    nfields = blast.setFormat('qid qlen qbegin qend sid slen sbegin send len pid evalue doc')
    for f in blast.fields:
        print('        ', f, '\t', getattr(blast, f))

    n=0
    while blast.next():
        n += 1
        print('   ', n, blast.line)
        print('   {}\t{}\t{}'.format(n, blast.sid, blast.evalue))

exit(0)
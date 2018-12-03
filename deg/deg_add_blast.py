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

    topmax = 3
    blastfile = sys.argv[1]
    sys.stderr.write('deg_add_blast')
    sys.stderr.write('\tBlast file: {}\n'.format(blastfile))

    blast = Blast()
    blast.new(blastfile)

    # format (diamond)
    nfields = blast.setFormat('qid qlen qbegin qend sid slen sbegin send len pid evalue doc')

    # TODO read in and store DEG file

    # read blast file, skipping any ID not in DEG file
    # TODO add logic for skipping


    n=0
    top = 1
    blast.read()
    query = blast.qid
    info = '\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(blast.qbegin,blast.qend,blast.sbegin,blast.send,blast.evalue,blast.doc)
    while blast.next():

        if blast.qid == query:
            top += 1

        else:
            sys.stderr.write('{}\t{}'.format(query, info))
            n += 1
            query = blast.qid

            top = 1
            info = ''

        if top <= topmax:
            info += '\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(blast.qbegin,blast.qend,blast.sbegin,blast.send,blast.evalue,blast.doc)

        if n > 5:
            break

exit(0)
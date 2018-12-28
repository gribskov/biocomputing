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
    degfile = sys.argv[2]
    outfile = 'merge.tsv'

    sys.stderr.write('deg_add_blast\n')
    sys.stderr.write('\tBlast file: {}\n'.format(blastfile))
    sys.stderr.write('\tDEG file: {}\n'.format(degfile))
    sys.stderr.write('\tOutput file: {}\n'.format(outfile))

    try:
        out = open(outfile, 'w')
    except:
        sys.stderr.write('Error opening output file ({})\n'.format(outfile))

    # open blastfile right away in case of error
    blast = Blast()
    blast.new(blastfile)
    # format (diamond)
    nfields = blast.setFormat('qid qlen qbegin qend sid slen sbegin send len pid evalue doc')

    # read in and store DEG file
    try:
        deg = open(degfile, 'r')
    except:
        sys.stderr.write('Unable to open DEG file ({})\n'.format(degfile))
        exit(1)

    ndeg = 0
    deginfo = {}
    header = deg.readline()
    for line in deg:
        # line = line.rstrip()
        ndeg += 1
        id, info = line.rstrip().replace('"', '').split('\t', maxsplit=1)
        deginfo[id] = info

    sys.stderr.write('\n{} genes read from {}\n'.format(ndeg, degfile))

    # read blast file, skipping any ID not in DEG file

    nblast = 0
    nfound = 0
    top = 1
    blast.read()
    query = blast.qid
    info = '{}\t{}\t{}\t{}\t{}\t{}'.format(blast.qbegin, blast.qend, blast.sbegin, blast.send,
                                           blast.evalue, blast.doc)
    while blast.next():

        if blast.qid == query:
            # another line of the current query
            top += 1

        else:
            # a new query, save the old query top lines
            if query in deginfo:
                nfound += 1
                out.write('{}\t{}\t"{}"\n'.format(query, deginfo[query], info))
            nblast += 1
            query = blast.qid
            info = '{}\t{}\t{}\t{}\t{}\t{}'.format(blast.qbegin, blast.qend, blast.sbegin,
                                                   blast.send,
                                                   blast.evalue, blast.doc)

            top = 1
            info = ''

        if top <= topmax:
            info += '\n{}\t{}\t{}\t{}\t{}\t{}'.format(blast.qbegin, blast.qend, blast.sbegin,
                                                      blast.send, blast.evalue, blast.doc)
        # if nfound > 5:
        #     break

    if query in deginfo:
        nfound += 1
        out.write('{}\t{}\t"{}"\n'.format(query, deginfo[query], info))

    sys.stderr.write('{} genes found in blast result {}\n'.format(nfound, blastfile))

exit(0)

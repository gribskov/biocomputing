"""-------------------------------------------------------------------------------------------------
deg_add_blast

Add blast search information to a tab delimited gene expression table.  The first column of the
deg table must match the query name of the blast search

usage

deg_add_blast.py <list file> <blast_search>

3 December 2018     Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys
from blast.blast import Blast


def match(deginfo, blast, level=None, filter="TRINITY_"):
    """-------------------------------------------------------------------------------------------------
    match the query ID with the names in the DEG file
    -------------------------------------------------------------------------------------------------"""
    qid = blast.qid.replace(filter, '')
    if level is None:
        # assume names match
        matchname = qid
    else:
        # assume a trinity name
        field = qid.split('_')
        matchname = '_'.join(field[:level + 1])
        blast.matchname = matchname

    if matchname in deginfo:
        return True
    return False


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
    except OSError:
        sys.stderr.write('Error opening output file ({})\n'.format(outfile))
        exit(1)

    # open blastfile right away in case of error
    blast = Blast()
    blast.new(blastfile)
    nfields = blast.setFormat(preset='diamond_doc')

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
    top = []
    matched_transcripts = []

    while blast.next():
        nblast += 1

        if match(deginfo, blast, level=1):
            # another line of the current query
            if len(top) == 0:
                # start a new list of hits
                nfound += 1
                query = blast.matchname
                matched_transcripts.append(query)

            if blast.matchname == query:
                # another match to the current query, add to list of top hits for this query
                top.append({'qid':    blast.qid, 'qlen': blast.qlen, 'qstart': blast.qstart, 'qend': blast.qend,
                            'sid':    blast.sid, 'slen': blast.slen, 'sstart': blast.sstart, 'send': blast.send,
                            'evalue': blast.evalue, 'doc': blast.doc.replace('UniRef50_', '').rstrip()})

            else:
                # a different query, write out the first topmax results
                top = sorted(top, key=lambda i: float(i['evalue']))[:topmax]
                out.write(f'{query}\t{deginfo[query]}\t')
                s = ''
                indent = '\t' * 12
                for n in range(len(top)):
                    t = top[n]
                    s += f"q:{t['qlen']},{t['qstart']},{t['qend']}\t"
                    s += f"s:{t['slen']},{t['sstart']},{t['send']}\t"
                    s += f"{t['evalue']} {t['doc']}"

                    n += 1
                    if n < len(top):
                        # s += f'\n{indent}'        # for multiple lines
                        s += '\t'

                out.write(f'{s}\n')
                top = []

    if top:
        nfound += 1
        top = sorted(top, key=lambda i: float(i['evalue']))[:topmax]
        out.write(f'{query}\t{deginfo[query]}\t')
        s = ''
        indent = '\t' * 12
        for n in range(len(top)):
            t = top[n]
            s += f"q:{t['qlen']},{t['qstart']},{t['qend']}\t"
            s += f"s:{t['slen']},{t['sstart']},{t['send']}\t"
            s += f"{t['evalue']} {t['doc']}"
            if n < len(top):
                # s += f'\n{indent}'        # for multiple lines
                s += '\t'

        out.write(f'{s}\n')

    # unmatched transcripts
    for id in deginfo:
        if id in matched_transcripts:
            continue

        out.write(f'{id}\t{deginfo[id]}\n')

    out.close()
    sys.stderr.write(f'{nblast} blast hits examined in {blastfile}\n')
    sys.stderr.write(f'{nfound} genes found')

    exit(0)

"""=================================================================================================
For a blast result comparing a set of trinity predicted transcripts to a database (probably
uniref50), Select the transcript with the greatest hit coverage at a specified level.
select level = 1 to 4, for a transcript named TRINITY_DN1371_c0_g2_i1
1 DN1371
2 DN1371_c0
3 DN1371_c0_g2
4 DN1371_c0_g2_i1

Michael Gribskov     13 November 2022
================================================================================================="""
import sys
from blast.blast import Blast


def read_genes(filename):
    """---------------------------------------------------------------------------------------------
    read in the list of genes
    :param filename: string
    :return: list of strings
    ---------------------------------------------------------------------------------------------"""
    select = open(filename, 'r')
    genes = []
    for line in select:
        name_in = line.rstrip()
        if line.find(' ') > -1:
            name_in = line.split(' ')[0]

        name_in.replace('TRINITY', '')
        genes.append(name_in)

    return genes


def match(deginfo, blast, level=None, filter="TRINITY_"):
    """-------------------------------------------------------------------------------------------------
    match the query ID with the names in the DEG file

    :param deginfo: list of string      names of selected queries
    :param blast: Blast object          search results
    :param level: int                   level for truncating trinity names
    :param filter: string               string to remove from blast query id (qid)
    :return: logical                    True if truncated query id is in deginfo
    -------------------------------------------------------------------------------------------------"""
    qid = blast.qid.replace(filter, '')
    if level is None:
        # assume names match
        matchname = qid
    else:
        # assume a trinity name
        field = qid.split('_')
        matchname = '_'.join(field[:level])
        blast.matchname = matchname

    if matchname in deginfo:
        return True
    return False


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    select_level = 2

    blastfile = sys.argv[1]  #
    degfile = sys.argv[2]  # list of transcript names , i.e. prefiltered set from DESeq2
    outfile = 'best_blast.tsv'

    sys.stderr.write('best_blast_from_trinity\n')
    sys.stderr.write('\tBlast file: {}\n'.format(blastfile))
    sys.stderr.write('\tSelected transcript file: {}\n'.format(degfile))
    sys.stderr.write('\tOutput file: {}\n\n'.format(outfile))

    try:
        out = open(outfile, 'w')
    except OSError:
        sys.stderr.write('Error opening output file ({})\n'.format(outfile))
        exit(1)

    # open blastfile right away in case of error
    blast = Blast()
    blast.new(blastfile)
    nfields = blast.setFormat(preset='diamond_doc')

    # read genes
    genes = read_genes(degfile)
    sys.stderr.write(f'\n{len(genes)} genes read\n')

    # find all matching blast hits, then write out the name of the best

    nblast = 0
    nmatch = 0
    nout = 0
    bestid = ''
    bestevalue = 100
    bestlen = 0
    bestsid = ''
    bestdesc = ''
    bestshort = ''
    query = ''
    sys.stderr.write(f'processing blast results in {blastfile}\n\n')
    while blast.next():
        nblast += 1
        sys.stderr.write(f'\r{nblast}')

        if match(genes, blast, level=select_level):
            nmatch += 1

            if blast.matchname == query:
                # matches to previous hit, update best

                isbest = False
                if float(blast.evalue) < bestevalue:
                    isbest = True
                elif float(blast.evalue) == bestevalue:
                    if int(blast.allen) > bestlen:
                        isbest = True

                if isbest:
                    bestid = blast.qid
                    bestevalue = float(blast.evalue)
                    bestsid = blast.sid
                    bestdesc = blast.doc
                    bestlen = int(blast.allen)
                    bestshort = blast.matchname


            # if blast.matchname == query:
            #     # another match to the current query, add to list of top hits for this query
            #     top.append({'qid':    blast.qid, 'qlen': blast.qlen, 'qstart': blast.qstart, 'qend': blast.qend,
            #                 'sid':    blast.sid, 'slen': blast.slen, 'sstart': blast.sstart, 'send': blast.send,
            #                 'evalue': blast.evalue, 'doc': blast.doc.replace('UniRef50_', '').rstrip()})

            else:
                # this is a new query, save the current result and reset best
                nout += 1
                if bestid:
                    out.write(f'{bestshort}\t{bestid}\t{bestsid}\t{bestevalue}\t{bestdesc}\n')
                bestid = blast.qid
                bestevalue = float(blast.evalue)
                bestlen = int(blast.allen)
                bestsid = blast.sid
                bestdesc = blast.doc
                bestshort = blast.matchname

            query = blast.matchname

    # end of loop over blast results
    out.write(f'{bestshort}\t{bestid}\t{bestsid}\t{bestevalue}\t{bestdesc}')  # don't forget the last one

    sys.stderr.write(f'\n\nBlast results examined: {nblast}\n')
    sys.stderr.write(f'ID matches: {nmatch}\n')
    sys.stderr.write(f'Best hits written: {nout}\n')

exit(0)

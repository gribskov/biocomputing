"""=================================================================================================
Create fasta sequences based on annotation. This could be used to get fasta sequences for stringtie
predicted genes

Michael Gribskov     10 March 2025
================================================================================================="""


def get_feature(gtf, target):
    """---------------------------------------------------------------------------------------------
    generator to extract the next matching feature from a gtf file
    TODO add GFF file
    format:
    ST4.03ch00	StringTie	transcript	5	2347	1000	+	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.1";

    :param gtf: fh          gtf file open for reading
    :param target: string  feature to select (column 2)
    :return: dict
    ---------------------------------------------------------------------------------------------"""
    feature_n = 0
    for line in gtf:
        if line.startswith('#'):
            continue

        sequence, source, feature, begin, end, score, strand, frame, annotation = line.rstrip().split('\t')
        if feature != target:
            continue

        begin = int(begin)
        end = int(end)
        yield {'sequence': sequence, 'source': source, 'feature': feature, 'begin': begin, 'end': end,
               'score':    score, 'strand': strand, 'frame': frame, 'annotation': annotation}

    return True

def get_fasta(fastafile):
    """---------------------------------------------------------------------------------------------
    Generator that returns the next fasta sequence from the fastafile.

    :param fastafile: string    path to fasta sequence file
    :return: dict               keys: id, documentation, sequence
    ---------------------------------------------------------------------------------------------"""
    fasta = open(fastafile, 'r')
    entry = {'id': '', 'documentation': '', 'sequence': ''}
    pos = 0
    for line in fasta:
        if line.startswith('>'):
            # id/documentation line
            if entry['sequence']:
                yield entry
                entry = entry = {'id': '', 'documentation': '', 'sequence': ''}
            field = line.rstrip().split(maxsplit=1)
            if len(field) > 1:
                entry['documentation'] = field[1]
            entry['id'] = field[0].replace('>','')

        else:
            # sequence line
            pos += len(line)
            print(pos)
            if pos>1000000:
                # TODO remove debug
                yield entry
            entry['sequence'] += line.rstrip()

    if entry['sequence']:
        yield entry

    fasta.close()
    return

def transcript_id(idstr):
    """---------------------------------------------------------------------------------------------
    extract transcript_id from a GTF annotation field, example
    gene_id "MSTRG.51021"; transcript_id "MSTRG.51021.1";
    :param idstr:
    :return:
    ---------------------------------------------------------------------------------------------"""
    field = idstr.split(';')
    for f in field:
        key, value = f.split(maxsplit=1)
        if key == 'transcript_id':
            return value.strip('"')
    return False

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gtffile = 'data/merged.gtf'
    fastafile = 'data/Stuberosum_448_v4.03.fa'
    transcript_outfile = 'data/transcript.out.fa'

    gtf = open(gtffile, 'r')
    feature = 'transcript'
    feature_n = 0
    transcripts = []
    for t in get_feature(gtf, feature):
        feature_n += 1
        transcripts.append(t)

    print(f'features read: {feature_n} from {gtffile}')
    gtf.close()

    out = open(transcript_outfile, 'w')
    for f in get_fasta(fastafile):
        seqid = f['id']
        # seqid='ST4.03ch01'
        for t in sorted( transcripts, key=lambda tr:(tr['sequence']==seqid, tr['begin'])):
            tid = transcript_id(t['annotation'])
            print(f"{tid}\t{t}")
            out.write(f">{tid} {seqid} begin:{t['begin']} end:{t['end']} strand:{t['strand']}\n")
            out.write(f"{f['sequence'][t['begin']-1:t['end']]}\n")

exit(0)

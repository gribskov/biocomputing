"""#################################################################################################
get annotation for transcripts from ensembl fa files
prototyped for mus-rel93

usage:
    anno_from_fasta.py <ref_base_name> <transcript_list>

15 August 2018     Michael Gribskov
#################################################################################################"""
import sys

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    base = sys.argv[1]
    transcriptfile = sys.argv[2]
    pepfile = base + ".pep.all.fa"

    sys.stderr.write('reference base name: {}\n'.format(sys.argv[1]))
    sys.stderr.write('peptide fasta file: {}\n'.format(pepfile))
    sys.stderr.write('transcript list: {}\n'.format(transcriptfile))

    # read annotations from peptide file

    try:
        pep = open(pepfile, 'r')
    except:
        sys.stderr.write('Unable to open petide fasta file ({})\n'.format(pepfile))

    try:
        transcript = open(transcriptfile, 'r')
    except:
        sys.stderr.write('Unable to open transcript fasta file ({})\n'.format(transcriptfile))

    # order of fields in annotation string
    order = ['pepid', 'gene', 'gene_biotype', 'gene_symbol',
             'release', 'chr', 'begin', 'end', 'dir']

    n = 0
    annotation = {}
    for line in pep:
        if line.startswith('>'):
            n += 1
            # sys.stderr.write(line)
            field = line.rstrip().split(' ')
            info = {}
            for i in range(len(field)):
                # sys.stderr.write('{} {}\n'.format(i, field[i]))
                if i == 1:
                    continue

                if i == 0:
                    # peptide id
                    info['pepid'] = field[i].lstrip('>')

                elif i == 2:
                    (tag,
                     info['release'],
                     info['chr'],
                     info['begin'],
                     info['end'],
                     info['dir']) = field[i].split(':', maxsplit=5)

                else:
                    (tag, value) = field[i].split(':', maxsplit=1)
                    if tag == 'description':
                        value = value + ' ' + ' '.join(field[i + 1:])
                        info['description'] = value
                        break

                    # all tags other than description
                    info[tag] = value

            if 'transcript' in info:

                anno_str = '{}\t'.format(info['transcript'])
                for k in order:
                    anno_str += '\t' + info[k]
                if 'description' in info:
                    anno_str += '\t' + info['description']

                annotation[info['transcript'].rstrip()] = anno_str

    pep.close()
    # done reading annotation

    # read gene list and match to transcripts

    for this_transcript in transcript:
        #this_transcript, isoform = t.split('.')
        this_transcript = this_transcript.rstrip()
        if this_transcript in annotation:
            sys.stdout.write('{}\t{}\n'.format(this_transcript, annotation[this_transcript]))
        else:
            sys.stdout.write('{} not found\n'.format(this_transcript))

    exit(0)

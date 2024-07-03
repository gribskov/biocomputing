"""=================================================================================================
get cDNA sequences from genome fa file based on Ensembl GTF

Michael Gribskov     03 July 2024
================================================================================================="""
import sys

sys.path.append('../gff')
from gtf2gff import Gtf

# translation function to complement bases
complement = str.maketrans('ACGTUacgtu', 'TGCAAtgcaa')


class Gene():
    """=============================================================================================
    Information about genes/transcripts

    ============================================================================================="""

    def __init__(self, gene_id='', seq='', begin=0, end=0, strand=''):
        self.gene_id = gene_id
        self.seq = seq
        self.begin = begin
        self.end = end
        self.strand = strand
        self.transcript = []

    def transcript_add(self, transcript_id='', begin=0, end=0):
        """-----------------------------------------------------------------------------------------
        add a new transcript to the list of transcripts for this gene

        :param transcript_id: string
        :param begin: int
        :param end: int
        :return: int                        number of transcripts in list
        -----------------------------------------------------------------------------------------"""
        self.transcript.append({'transcript_id': transcript_id,
                                'begin':         begin,
                                'end':           end,
                                'exon':          []
                                })
        return len(self.transcript)

    def exon_add(self, begin=0, end=0):
        """-----------------------------------------------------------------------------------------
        Add an exon to the last transcript in the list
        :param begin:
        :param end:
        :return:
        -----------------------------------------------------------------------------------------"""
        self.transcript[-1]['exon'].append({'begin': begin, 'end': end})
        return len(self.transcript[-1])


def get_spliced(target, seq):
    """---------------------------------------------------------------------------------------------
    Extract the spliced sequence for each transcript in t from the current sequence
    
    :param target: Gene         transcript sequences to extract
    :param seq: string          genomic sequence
    :return: list of string     constructed spliced sequences
    ---------------------------------------------------------------------------------------------"""
    spliced = []
    for t in target.transcript:
        ss = ''
        for e in t['exon']:
            ss += seq[e[0] - 1:e[1]]
        if t.strand == '-':
            # reverse complement
            ss = complement(ss[:-1])

        spliced.append(ss)

    return spliced


def write_out(out, target, spliced, linelen=100):
    """---------------------------------------------------------------------------------------------
    Write the spliced sequences to the output file

    :param out: file handle             file open for writing
    :param target: Gene                 transcript sequences to extract
    :param spliced: list of string      constructed spliced sequences
    :return: True
    ---------------------------------------------------------------------------------------------"""
    n = 0
    for t in target.transcript:
        out.write(f">{t['transcript_id']} gene:{target.gene_id} {target.seq}:{t['begin']}:{t['end']}\n")
        pos = 0
        while pos < len(spliced[n]):
            out.write(f'{spliced[n][pos:pos + linelen]}\n')
        pos += linelen
        n += 1

    return


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    genomefile = 'A:/mrg\Dropbox/21dog_dlbcl/reference data/cfam112/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa'
    try:
        genome = open(genomefile, 'r')
    except OSError:
        sys.stderr.write(f'Unable to open genome file ({genomefile}')
        exit(1)

    gtf_file = 'A:/mrg\Dropbox/21dog_dlbcl/reference data/cfam112/Canis_lupus_familiaris.ROS_Cfam_1.0.112.gtf'
    gtf = Gtf(gtf_file)

    transcript_list = []
    nt = 0
    while gtf.next():
        nt += 1
        if nt > 50000:
            break
        gtf.parse()
        gtf.add_ensemble_attributes()

        parsed = gtf.parsed
        if parsed['feature'] == 'gene':
            print(f"{parsed['seqname']}\t{parsed['start']}\t{parsed['end']}\t{parsed['strand']}\t{parsed['gene_id']}")
            gene = Gene(parsed['gene_id'], parsed['seqname'], parsed['start'], parsed['end'])
            transcript_list.append(gene)

        elif parsed['feature'] == 'transcript':
            print(f"\t{parsed['start']}\t{parsed['end']}\t{parsed['transcript_id']}")
            gene.transcript_add(parsed['transcript_id'], parsed['start'], parsed['end'])

        # utrs are included in exons
        # elif parsed['feature'] in ('exon', 'three_prime_utr', 'five_prime_utr'):
        elif parsed['feature'] in ('exon'):
            print(f"\t\t{parsed['feature']}\t{parsed['start']}\t{parsed['end']}")
            gene.exon_add(parsed['start'], parsed['end'])

    print(f'Counting data')
    transcript_n = 0
    exon_n = 0
    chr_idx = {}
    for t in sorted(transcript_list, key=lambda x: (x.seq, x.begin)):
        if t.seq in chr_idx:
            chr_idx[t.seq].append(t)
        else:
            chr_idx[t.seq] = [t]

        print(f'chromosome {t.seq}\t {t.begin}')
        transcript_n += len(t.transcript)
        for this_transcript in t.transcript:
            exon_n += len(this_transcript['exon'])

    print(f'{len(chr_idx)} sequences found')
    print(f'{transcript_n} transcripts')
    print(f'{exon_n} exons')

    # read the sequence file one sequence one line at a time. Once we know the chromosome, start looking for a chunk
    # with the next transcript (in sorted order in chr_idx). this way only a fraction of the sequence is in memory

    out = open('transcripts.fa', 'w')

    current_id = ''
    current_doc = ''
    seq = ''
    nline = 0
    for line in genome:
        print(f'{nline}:{len(seq)}\n{line}')
        nline += 1

        if line.startswith('>'):
            # title line of new sequence
            id, doc = line.rstrip().split(' ', maxsplit=1)
            id = id.replace('>', '')
            if current_id != id and seq:
                # new id and sequence exists, process all ranges
                # content of chr_idx are  Gene objects
                for t in chr_idx[current_id]:
                    print(f'{t.seq}:{t.gene_id}\t{t.begin}\t{t.end}')
                    spliced = get_spliced(t, seq)
                    write_out(out, t, spliced)

            if current_id != id:
                # new sequence
                seq = ''
                current_id = id
                current_doc = doc

        else:
            seq += line.rstrip()

    # don't forget the last sequence

    out.close()
    exit(0)

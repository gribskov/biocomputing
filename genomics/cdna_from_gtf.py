"""=================================================================================================
get cDNA sequences from genome fa file based on Ensembl GTF

Michael Gribskov     03 July 2024
================================================================================================="""
import sys

sys.path.append('../gff')
from gtf2gff import Gtf


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


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    gtf_file = 'A:/mrg\Dropbox/21dog_dlbcl/reference data/cfam112/Canis_lupus_familiaris.ROS_Cfam_1.0.112.gtf'
    gtf = Gtf(gtf_file)

    transcript_list = []
    while gtf.next():
        gtf.parse()
        gtf.add_ensemble_attributes()

        parsed = gtf.parsed
        if parsed['feature'] == 'gene':
            print(f"{parsed['seqname']}\t{parsed['start']}\t{parsed['end']}\t{parsed['strand']}\t{parsed['gene_id']}")
            gene = Gene(parsed['seqname'], parsed['gene_id'], parsed['start'], parsed['end'])
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
    sequences = []
    transcript_n = 0
    exon_n = 0
    for t in transcript_list:
        if t.seq not in sequences:
            sequences.append( t.seq)
        transcript_n += len(t.transcript)
        for this_transcript in t.transcript:
            exon_n += len(this_transcript['exon'])

    print(f'{len(sequences)} sequences found')
    print(f'{transcript_n} transcripts')
    print(f'{exon_n} exons')

    exit(0)

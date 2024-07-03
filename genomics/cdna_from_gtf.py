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


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    gtf_file = 'A:/mrg\Dropbox/21dog_dlbcl/reference data/cfam112/Canis_lupus_familiaris.ROS_Cfam_1.0.112.gtf'
    gtf = Gtf(gtf_file)

    transcript_list = []
    while gtf.next():
        gtf.parse()
        # gtf.attribute2gff()
        gtf.add_ensemble_attributes()

        parsed = gtf.parsed
        if parsed['feature'] == 'gene':
            print(f"{parsed['seqname']}\t{parsed['start']}\t{parsed['end']}\t{parsed['strand']}\t{parsed['gene_id']}")
            gene = Gene(parsed['seqname'], parsed['gene_id'], parsed['start'], parsed['end'])
            transcript_list.append(gene)
        elif parsed['feature'] == 'transcript':
            print(f"\t{parsed['start']}\t{parsed['end']}\t{parsed['transcript_id']}")
            gene.transcript.append({'transcript_id': parsed['transcript_id'],
                                    'begin':         parsed['start'],
                                    'end':           parsed['end'],
                                    'exon':          []
                                    })
            exon = gene.transcript[-1]['exon']
        elif parsed['feature'] in ('exon', 'three_prime_utr', 'five_prime_utr'):
            print(f"\t\t{parsed['feature']}\t{parsed['start']}\t{parsed['end']}")
            exon.append({'begin': parsed['start'], 'end': parsed['end']})

    exit(0)

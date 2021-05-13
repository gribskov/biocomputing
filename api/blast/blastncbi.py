"""=================================================================================================
Run Blast queries at the NCBI server

Michael Gribskov     31 March 2021
================================================================================================="""
import sys
import time
import requests
from bs4 import BeautifulSoup, Comment
from ..jobmanager_api import JobManagerAPI


class BlastNCBI(JobManagerAPI):
    """=============================================================================================
    Not all options are implemented.  See https://ncbi.github.io/blast-cloud/dev/api.html
    Some options are taken directly from the web form

    ============================================================================================="""
    available = {'program':  ['blastp', 'psiBlast', 'deltaBlast', 'kmerBlastp', 'phiBlast',
                              'blastx',
                              'blastn', 'megablast', 'discomegablast', 'tblastn', 'tblastx'],

                 'database': ['nr', 'refseq_select_prot', 'refseq_protein', 'SMARTBLAST/landmark',
                              'swissprot', 'pataa', 'pdb', 'env_nr', 'tsa_nr', 'nt',
                              'refseq_select_rna', 'refseq_rna', 'refseq_representative_genomes',
                              'refseq_genomes', 'Whole_Genome_Shotgun_contigs', 'est', 'sra',
                              'htgs', 'patnt', 'pdbnt', 'genomic/9606/RefSeqGene', 'gss', 'sts',
                              'rRNA_typestrains/16S_ribosomal_RNA', 'TL/18S_fungal_sequences'
                                                                    'TL/28S_fungal_sequences',
                              'rRNA_typestrains/ITS_RefSeq_Fungi'],

                 'matrix':   ['PAM30', 'PAM70', 'PAM250', 'BLOSUM80', 'BLOSUM62', 'BLOSUM45',
                              'BLOSUM50', 'BLOSUM90'],

                 'format':   ['HTML', 'Text', 'XML', 'XML2', 'JSON2', 'Tabular']
                 }

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""

        self.url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
        self.program = ''
        self.database = ''
        self.entrez = ''
        self.query = ''
        self.email = ''
        self.format = 'XML'
        self.expect = None  # expect - E-value threshold
        self.num_seq = None  # NUM_SEQ - number of hits to report
        self.matrix = ''  # matrixName - scoring matrix
        self.qorganism = 1  # qorganism - TaxID to limit query
        self.entrez = ''  # qquery

        self.jobstatus = ''
        self.jobid = ''  # NCBI calls this RID
        self.rtoe = None  # estimated time for job
        self.response = None

        # jobname, message, and content inherited from JobManagerAPI

    def get_qblastinfo(self):
        """-----------------------------------------------------------------------------------------
        The blast server returns status information in a comment labelled QBlastInfo.  Returns a
        dictionary with the tags and fields present in this block

        RID = request ID
        RTOE = Request Time Of Execution = estimated run time
        Status - for a running job, one of WAITING, UNKNOWN, or READY

        :return: dict, tags and alues inside the QBlastInfo block (RID and RTOE, or Status)
        -----------------------------------------------------------------------------------------"""
        # To parse correctly, there must be whitespace between the <!-- and 'QBlastInfoBegin'
        soup = BeautifulSoup(self.response.text.replace('QBlastInfoBegin', ' QBlastInfoBegin'),
                             'html.parser')
        comment = soup.find(
            text=lambda text: isinstance(text, Comment) and text.find('QBlastInfoBegin') > 0)

        info = {}
        # field = comment.split('\n')
        for line in comment.split('\n'):
            line = line.strip()
            if not line or line.startswith('QBlast'):
                continue

            tag, value = line.split('=')
            info[tag.strip()] = value.strip()

        return info

    def parse_xml(self):
        """-----------------------------------------------------------------------------------------
        Parse the xml response and return as a list of dictionaries.  Each hit may have multiple
        HSPs.

        :return: list of dictionaries
        -----------------------------------------------------------------------------------------"""
        soup = BeautifulSoup(self.response.text, "xml")
        # print(soup.prettify())

        hits = soup.find_all('Hit')

        blasthits = []
        for hit in hits:

            this_hit = {'id':        hit.Hit_id.text, 'def': hit.Hit_def.text,
                        'accession': hit.Hit_accession.text, 'len': hit.Hit_len.text,
                        'hsp':       []}
            blasthits.append(this_hit)
            for hsp in hit.find_all('Hsp'):
                hsp_num = hsp.find('Hsp_num').text
                this_hsp = {'hsp_num': hsp_num,
                            'bitscore':  hsp.find('Hsp_bit-score').text,
                            'hspscore':  hsp.find('Hsp_score').text,
                            'evalue':    hsp.find('Hsp_evalue').text,
                            'identity':  hsp.find('Hsp_identity').text,
                            'positive':  hsp.find('Hsp_positive').text,
                            'gaps':      hsp.find('Hsp_gaps').text,
                            'align_len': hsp.find('Hsp_align-len').text,
                            'qbegin':    hsp.find('Hsp_query-from').text,
                            'qend':      hsp.find('Hsp_query-to').text,
                            'sbegin':    hsp.find('Hsp_hit-from').text,
                            'send':      hsp.find('Hsp_hit-to').text}
                this_hit['hsp'].append(this_hsp)

        return blasthits

    def result(self):
        """-----------------------------------------------------------------------------------------
        Retrieve the result from the server

        :return: string, status of retrieval
        -----------------------------------------------------------------------------------------"""

        command = {'CMD':         'Get',
                   'FORMAT_TYPE': self.format,
                   'RID':         self.jobid
                   }
        self.response = requests.post(self.url, command)
        self.content = self.response.text
        self.message = {'type':     'retrieved',
                        'text':     f'job_id={self.jobid};output_len={len(self.content)}',
                        'loglevel': 2}
        # TODO add check for failure

        return 'retrieved'

    def status(self):
        """-----------------------------------------------------------------------------------------
        Polls the server and gets the status of the job in self.jobid. The server returns:
        WAITING, UNKNOWN, or READY
        self.status is set as running, unknown, or finished, respectively

        :return: string, status running | unknown | finished
        -----------------------------------------------------------------------------------------"""
        xlate = {'WAITING': 'running', 'UNKNOWN': 'unknown', 'READY': 'finished'}

        command = {'CMD':           'Get',
                   'FORMAT_OBJECT': 'SearchInfo',
                   'RID':           self.jobid
                   }
        self.response = requests.post(self.url, command)
        # print(self.response.request.body, '\n')

        info = self.get_qblastinfo()
        self.jobstatus = xlate[info['Status']]
        self.message = {'type':     self.jobstatus,
                        'text':     f'job_id={self.jobid};response={info["Status"]}',
                        'loglevel': 2}

        return self.jobstatus

    def submit(self):
        """-----------------------------------------------------------------------------------------
        submit query to server at NCBI

        :return: string, request ID (rid)
        -----------------------------------------------------------------------------------------"""
        command = {'CMD':      'Put',
                   'PROGRAM':  self.program,
                   'DATABASE': self.database,
                   'QUERY':    self.query,
                   'EMAIL':    self.email
                   }
        # print('command:', command)
        self.response = requests.post(self.url, command)
        info = self.get_qblastinfo()

        self.jobid = info['RID']
        self.rtoe = int(info['RTOE']) * 60
        self.jobstatus = 'submitted'
        self.message = {'type':     'submitted',
                        'text':     f'job_name={self.jobname};job_id={self.jobid}',
                        'loglevel': 1}

        return self.jobid


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    blast = BlastNCBI()
    blast.email = 'gribskov@purdue.edu'
    blast.program = 'blastp'
    blast.database = 'pdb'
    #     blast.query = '''>HA
    # MGEIGFTEKQEALVKESWEILKQDIPKYSLHFFSQILEIAPAAKGLFSFLRDSDEVPHNNPKLKAHAVKV
    # FKMTCETAIQLREEGKVVVADTTLQYLGSIHLKSGVIDPHFEVVKEALLRTLKEGLGEKYNEEVEGAWSQ
    # AYDHLALAIKTEMKQEES
    # '''
    blast.query = '''>HA
MGEIGFTEKQEALVKESWEILKQDIPKYSLHFFSQILEIAPAAKGLFSFLRDSDEVPHNNPKLKAHAVKV
SFDARPVDAQAAVDDAADADAPPPASDGALVNERTPLLHTASHHSAEGDVTPPSEEWSSHALLRTLKEG
GRDPEMGEKKILELVQACDEWLELPPRDLEKPFLMPVEDVFSISGRGTVATGRVERGIAT
THRLESKSQVNLIDTLRKQLEEAGPLINAVRASSALMETDAAKQKMENEKLQAELDKAKA
LGEKYNEEVEGAWSQAYDHLALAIKTEMKQEES
'''
    jobid = blast.submit()
    print('rid:', jobid, '\trtoe:', blast.rtoe)
    poll_delay = 60
    while blast.jobstatus != 'finished':
        time.sleep(poll_delay)
        event_time = time.strftime('%d/%b/%G:%H:%M:%S', time.localtime(time.time()))
        print(event_time)
        blast.status()
        print(blast.jobstatus)

    blast.result()
    print(blast.response.text)

    exit(0)

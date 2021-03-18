import sys
import time
import json
import requests

class Interpro:
    """=============================================================================================
    Interpro class for running interproscan

    25 December 2018    Michael Gribskov
    ============================================================================================="""

    def __init__(self, loglevel=0, poll_time=60, poll_count=20):
        """-----------------------------------------------------------------------------------------
        interpro query/response constructor

        loglevel   0 no log, 1 job submission/completion, 2 all
        -----------------------------------------------------------------------------------------"""
        self.log = loglevel
        self.log_fh = sys.stderr

        # availble options/parameters taken from
        # https://www.ebi.ac.uk/Tools/services/rest/iprscan5/parameterdetails/appl
        self.applications_avail = ['TIGRFAM', 'SFLD', 'Phobius', 'SignalP', 'SignalP_EUK',
                                   'SignalP_GRAM_POSITIVE', 'SignalP_GRAM_NEGATIVE', 'SUPERFAMILY',
                                   'Panther', 'Gene3d', 'HAMAP', 'PrositeProfiles',
                                   'PrositePatterns', 'Coils', 'SMART', 'CDD', 'PRINTS', 'PfamA',
                                   'MobiDBLite', 'PIRSF', 'TMHMM', ]

        self.commands_avail = ['run', 'status', 'result']

        # available outputs
        # from https://www.ebi.ac.uk/Tools/services/rest/iprscan5/resulttypes
        # /iprscan5-R20210318-175621-0074-61690295-p2m
        # log - The output from the tool itself
        # out - The results of the job (XML format)
        # tsv - The results of the job in text format, tab separated values
        # xml - The results of the job in XML
        # gff - The results of the job in GFF3 format
        # json - The results of the job in JSON format
        # htmltarball - The results of the job in a tarball zip file
        # sequence - Input sequence as seen by the tool
        # submission - The submission details which was submitted as a job
        self.output_avail = {'out', 'log', 'tsv', 'xml', 'gff', 'json',
                             'htmltarball', 'sequence', 'submission'}
        self.poll_time = poll_time  # seconds between polling
        self.poll_count = poll_count  # maximum number of times to poll

        self.email = ''  # user email (optional)
        self.title = ''  # title for job (optional)
        self.sequence = ''
        self.applications = []
        self.output = ''
        self.parameters = {}

        self.url = u'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/'
        self.jobid = ''
        self.state = 'UNKNOWN'
        self.content = ''

    def application_select(self, selected, keep=False):
        """-----------------------------------------------------------------------------------------
        Add a list of applications to be run.  Each application in the list is compared to the
        available applications and if not present a warning is issued.  The default is to run all
        applications so an empty
        list signifies the default.

        :param selected: list of strings, selected applications to run
        :param keep: Boolean, retain current applications, just add new ones
        :return: int, number of selected applications
        -----------------------------------------------------------------------------------------"""
        if not keep:
            self.applications = []

        for app in selected:
            if app in self.applications_avail:
                self.applications.append(app)
            else:
                self.log_fh.write(
                    '{}\trequested application ({}) not available'.format(Interpro.logtime(), app))

        return len(self.applications)

    def output_select(self, selected):
        """-----------------------------------------------------------------------------------------
        select the output format.  Only one is allowed

        :param selected: string, one of the formats in self.output_avail
        :return: True if format is available
        -----------------------------------------------------------------------------------------"""
        self.output = ''
        if selected in self.output_avail:
            self.output = selected
        else:
            self.log_fh.write(
                '{}\trequested output format ({}) is not available'.format(
                    Interpro.logtime(), selected))
            return False

        return True

    def parameter_select(self, select):
        """-----------------------------------------------------------------------------------------
        Select additional tag value pairs to add to parameters.  There is no checking so be correct

        :param select: dict
        :return: int number of parameters in dictionary
        -----------------------------------------------------------------------------------------"""
        for key in select:
            self.parameters[key] = select[key]

        return len(self.parameters)

    def run(self):
        """-----------------------------------------------------------------------------------------
        Construct a REST command and dispatch the job to the server
        Any previously existing jobID is overwritten
        :return: logical, True = success, False = failure
        -----------------------------------------------------------------------------------------"""
        is_success = False

        # general fields for all queries
        param = {u'email': self.email, u'title': self.title, u'sequence': self.sequence,
                 u'output':self.output}

        if self.applications:
            # add selected applications
            param['appl'] = ','.join(self.applications)


        if self.parameters:
            for para in self.parameters:
                param[para] = self.parameters[para]

        command = self.url + 'run'
        response = requests.post(command, files=param, headers={'User-Agent':'ips-client'})

        print(response.request.headers,'\n')
        print(response.request.body,'\n')

        if not self.response_is_error('submitting job', response):
            # success
            self.jobid = response.text
            if self.log:
                # TODO add sequence name?
                self.log_fh.write(
                    '{}\tinterproscan job {} submitted\n'.format(Interpro.logtime(), self.jobid))
            is_success = True

        return is_success


    def status(self):
        """-----------------------------------------------------------------------------------------
        Poll job status at the server.
        wait for self.poll_time seconds after polling

        :return: Logical, True if a result was returned
        -----------------------------------------------------------------------------------------"""
        response = None
        complete = False
        tries = 0
        while not complete:
            tries += 1

            command = self.url + 'status/' + self.jobid
            response = requests.get(command)
            if self.log > 1:
                self.log_fh.write('{}\tinterproscan job {} polling - {}\n'.format(
                    Interpro.logtime(), self.jobid, response.text))

            if 'FINISHED' in response.text:
                complete = True
                break
            elif tries >= self.poll_count:
                break

            # don't poll too often
            time.sleep(self.poll_time)

        self.state = response.text
        if not complete:
            # polling reached limit
            if self.log > 0:
                self.log_fh.write('{}\tinterproscan job {} {} - poll_limit={}\n'.format(
                    Interpro.logtime(), self.jobid, self.state, self.poll_count))

        else:
            if self.log > 0:
                self.log_fh.write(
                    '{}\tinterproscan job {} finished\n'.format(Interpro.logtime(), self.jobid))

        return complete


    def result(self):
        """-----------------------------------------------------------------------------------------
        Retrieve the result

        :return: Logical True=success, False=failure
        -----------------------------------------------------------------------------------------"""
        # get the final result
        command = self.url + 'result/' + self.jobid + '/' + self.output
        response = requests.get(command)
        if not self.response_is_error('retrieving result', response):
            # success
            self.content = response.text
            if self.log > 1:
                self.log_fh.write('{}\tinterproscan job {} result retrieved from {} as {}'.format(
                    Interpro.logtime(), self.jobid, self.url, self.output))
            return True

        return False


    def response_is_error(self, task, response):
        """-----------------------------------------------------------------------------------------
        Return true if the response code is other than 200. Write error message to stderr if
        loglevel > 1. Task is a string describing the task that failed for inclusion in the error
        message

        :param response: requests object response
        :return: logical True = error, False = no error
        -----------------------------------------------------------------------------------------"""
        is_error = False
        if not response.status_code == 200:
            if self.log > 0:
                self.log_fh.write('{}\t{}interproscan job {} error\tstatus={}'.format(
                    Interpro.logtime(), self.jobid, task, response.status_code))

            is_error = True

        return is_error


    def set_log_fh(self, fh):
        """-----------------------------------------------------------------------------------------
        The output for the log is STDERR by default.  This function allows you to change it
        -----------------------------------------------------------------------------------------"""
        self.log_fh = fh
        return fh


    @classmethod
    def logtime(cls):
        """-----------------------------------------------------------------------------------------
        Return current time as a string. Format is 10/Oct/2000:13:55:36 which is similar to the
        common log format (without the time zone)

        :return: string
        -----------------------------------------------------------------------------------------"""
        return time.strftime('%d/%b/%G:%H:%M:%S', time.localtime(time.time()))


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    protein_test = '''>AAX90616.1 Src [Mus musculus]
MGSNKSKPKDASQRRRSLEPSENVHGAGGAFPASQTPSKPASADGHRGPSAAFVPPAAEPKLFGGFNSSD
TVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYI
PSNYVAPSDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHY
KIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVCPTSKPQTQGLAKDAWEIPRESLRLEVK
LGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMN
KGSLLDFLKGETGKYLRLPQLVDMSAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIE
DNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMP
CPPECPESLHDLMCQCWRKEPEERPTFEYLQAFLEDYFTSTEPQYQPGENL
'''
    testpro='''>QAB35645.1 BRI1 protein [Solanum tuberosum]
MKAHKTVFYQHPLSLNKLFFVLLLIFFLPPASPASVNGLFKDSQQLLSFKAALPPTPTLLQNWLPSTDPC
SFTGVSCKNSRVSSIDLSNTFLSVDFSLVTSYLLPLSNLESLVLKNANLSGSLTSAAKSQCGVSLDSIDL
AENTISGPISDISSFGVCSNLKSLNLSKNFLDPPGKEILKGATFSLQVLDLSYNNISGFNLFPWVSSMGF
GELEFFSLKGNKLAGSIPELDFKNLSHLDLSANNFSTVFPSFKDCSNLQHLDLSSNKFYGDIGSSLSSCG
KLSFLNLTNNQFVGLVPKLQSESLQYLYLRGNDFQGVYPNQLADLCKTVVELDLSYNNFSGMVPESLGEC
SSLELVDISNNNFSGKLPVDTLLKLSNMKTMVLSFNKFVGVLPDSFSNLLKLETLDVSSNNLTGVIPSGI
CKDPMNNLKVLYLQNNLFEGPIPDSLSNCSQLVSLDLSFNYLTRRIPSSLGSLSKLKDLILWLYQLSGEI
PQELMYLQALENLILDFNDLTGPIPASLSNCTKLNWISLSNNQLSGEIPASLGRLSNLAILKLGNNSISG
NIPAELGNCQSLIWLDLNTNFLSGSIPPPLFKQSGNIAVALLTGKRYVYIKNDGSKECHGAGNLLEFGGI
RQEQLGRISTRHPCNFTRVYRGITQPTFNHNGSMIFLDLSYNKLEGSIPKELGTMYYLSILNLGHNDLSG
MIPQDHGGLKNVAILDLSYNRFNGPIPNSLTSLTLLGEIDLSNNNLSGMIPESAPFDTFPDYRFANNSLC
GYPLPLPCSSGPKSDANQHQKSHRRQASLAGSVAMGLLFSLFCIFGLIIVAIETKKRRKKKEAALEAYMD
GHSHSATANSAWKFTSAREALSINLAAFEKPLRKLTFADLLEATNGFHNDSLVGSGGFGDVYKAQLKDGS
VVAIKKLIHVSGQGDREFTAEMETIGKIKHRNLVPLLGYCKVGEERLLVYEYMKYGSLEDVLHDRKKIGI
KLNWPARRKIAIGAARRLAFLHHNCIPHIIHRDMKSSNVLLDENLEARVSDFGMARLMSAMDTHLSVSTL
AGTPGYVPPEYYQSFRCSTKGDVYSYGVVLLELLTGKQPTDSADFGDNNLVGWVKLHAKGKITDVFDREL
LKEDPSIEIELLQHLKVACACLDDRHWKRPTMIQVMAMFKEIQAGSGMDSTSTIGADDVNFSAVEGGIEM
GINESIKEGNELSKHL'''

    nucleotide_test = '''>AB512659.1 Bos javanicus birmanicus x Bos indicus HBB gene for hemoglobin beta, complete cds, allele: HBB-B, specimen_voucher: personal:Cambodia Hybrid-CamII-2
ACAAACAGACACCATGCTGACTGCTGAGGAGAAGGCTGCCGTCACCGCCTTTTGGAGCAAGGTGCATGTG
GATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAGGTATCCCACTTACAAGGCAGGTTTAAGGAGAGTGAAA
TGCACCTGGGCGTGTGAGGACAGAGCCGTCCCTGAGATTCTGAAAGCTGCTGGCTTCCTCTGACCTTGTG
CTGTTTTCTCCCCCTAGGCTGCTGGTTGTCTACCCCTGGACTCAGAGGTTCTTTGAGTCCTTTGGGGACT
TGTCCACTGCTGATGCTGTTATGAACAACCCTAAGGTGAAGGCCCATGGCAAGAAGGTGCTAGATTCCTT
TAGTAATGGCATGAAGCATCTCGATGACCTCAAGGGCACCTTTGCTGCGCTGAGTGAGCTGCACTGTGAT
AAGCTGCATGTGGATCCTGAGAACTTCAAGGTGAGTTTGTGGAATCCTCAGTGTTCTTCTTCTTTTTATG
GTCAAGCTCATGTCATGAGGAGAAAGCTGAATGGCAGGACACAGTTTAGAATGGAGAAGAGGTATTCTGG
TTAGATTACTAAGGACTCCTCAGAACCGTTTAGACTCTTTTAACCTCTTTGTTCACAACCAGCATTTCCT
CTGATTCATTCTTGTTCTCTGTTGTCTGCAATGTCCTCTTTTTAATTATATTTTTTATTTTGAGGGTTTA
ATTTGAAAAAAAATTATATATCAACTTTAAAAATTGTATCTAATATTTCCCCCTTATCTGTTTCTTTCAA
GGAATAAAATGTTCTATTGCTTTTTGAAATGATTCAAAATAATAAAAATAATAACAAGTTCTGGATTAAG
TTAGAAAGAGAGAAACATTTCTAAATATATATTCAGGAAGATATAGGTAGATTCACATCAGTAGTAACAA
CTTCACTTCAGTCATCTTTGTGCTTATATCCTACGGTCACAGCTTGGGATAAGACTGAAATACCCTGAAT
CTAACCTTGGATTTCCCTCATAGCTCAGTTGGTTAAGCATCTGCCTGCAATGCAAGAGATCCCAGTTCGA
TTCCTGGGTCGGGAAGGATGGCTGGAGAAGGGATAGGCACCCACTCCAGTATTCTTGGGTTTCCCTTGTG
GCTCAGCTGGTAAAGAATCTGCCTGCAATGTGGGAGACCCAGCTTCTATCCCTGAGTTTGGAGGATCCCC
TGGAGAAGGGAAAGGCTACCCACTCCAGTATTCTGGCCTGGAGAAATCTATGGACTGTATAGTCCATGGG
GTTGCAAAGAATCAGACACGATTGAGAGACTCTCACTTCACTCACCTGCACTAACCCTGCCCTTGCTTAA
TGTCTTTTCCACACAGCTCCTGGGCAACGTGCTAGTGGTTGTGCTGGCTCGCAATTTTGGCAACGAATTC
ACCCCGGTGCTGCAGGCTGACTTTCAGAAGGTGGTGGCTGGTGTGGCCAATGCCCTGGCCCACAGATATC
ATTAAGCTCCCTTTCCTGCTTTCCAGGAAAGGTTTTTTCATCCTCAGAGCCCAAAAATTGAATATGGAAA
AATTATGAAGTGTTTTGAGCATCTGGCCTCTGCCTAATAAAGACATTTATTTTCATTGCACTGGTGTATT
TAAATTATTTCACTGTCTCTTACTCAGATGGGCACATGGGAGGGCAAAACACTGAAGACATAAAGAAATG
AAGG
'''
    ips = Interpro(loglevel=2, poll_time=60)
    ips.email = 'gribskov@purdue.edu'
    ips.title = 'globin'
    ips.sequence = testpro
    # ips.application_select(['TIGRFAM', 'CDD', 'PfamA'])
    # ips.application_select(['Phobius', 'ProSitePatterns'])
    ips.output_select('json')
    ips.parameter_select({'goterms':True, 'pathways':True})

    if not ips.run():
        exit(1)

    if ips.status():
        ips.result()
        # sys.stdout.write(ips.content)
        j = json.loads(ips.content)

    print('done')

    exit(0)

import sys
import time
import json
import requests


class Interpro:
    """=============================================================================================
    Interpro class for running interproscan

    25 December 2018    Michael Gribskov
    ============================================================================================="""

    def __init__(self, loglevel=0, poll_time=60, poll_max=20):
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
        self.poll_max = poll_max  # maximum number of times to poll
        self.poll_count = 0  # number of times this job has been polled

        self.email = ''  # user email (optional)
        self.title = ''  # title for job (optional)
        self.sequence = ''
        self.applications = []
        self.output = 'json'
        self.parameters = {}

        self.url = u'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/'
        self.jobid = ''
        self.jobstatus = ''

        self.response = None
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
            if app == 'Pfam':
                app = 'PfamA'
            if app in self.applications_avail:
                self.applications.append(app)
            else:
                self.log_message('not_available', 'application={}'.format(app))

        return len(self.applications)

    def output_select(self, selected):
        """-----------------------------------------------------------------------------------------
        select the output format.  Only one is allowed

        :param selected: string, one of the formats in self.output_avail
        :return: True if format is available
        -----------------------------------------------------------------------------------------"""
        # self.output = ''
        if selected in self.output_avail:
            self.output = selected
        else:
            self.log_message('not_available', 'output={}'.format(selected))

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

    def parse_json(self):
        """-----------------------------------------------------------------------------------------
        Parse the contented returned from the server in JSON format. Three outputs are produced
            A list of dictionaries for each hit in the sequence
            A list of dictionaries listing  GO terms and what entries they were drawn from
                keys: 'name': gene ontology ID
                      'category': ontology category = BIOLOGICAL_pROCESS |
                                                      MOLECULAR_FUNCTION |
                                                      CELLULAR_COMPONENT,
                      'source': list of strings, motifs that are associated with this term
            A list of dictionaries listing pathways and the entry theies were drawn from
                keys: 'databaseName': the pathway database KEGG | Reactome | Metacyc
                      'id': UID of the pathway in its database
                      'name': text description of pathway
                      'source': list of strings, motifs that are associated with this term

        This is fairly specific for my purpose

        :return:
        -----------------------------------------------------------------------------------------"""
        pjson = json.loads(self.content)

        matches = pjson['results'][0]['matches']
        go_all = {}
        path_all = {}
        motifs = []
        for m in matches:
            # each match in the interproscan search
            # source_accession is the UID of the matching motif in the source database
            # source is the database and version listed at interpro
            # entry is the subtree of information about a match
            entry = m['signature']['entry']
            signature = m['signature']['signatureLibraryRelease']
            source_accession = m['signature']['accession']
            source = '{} {}'.format(signature['library'],
                                    signature['version'])

            if not entry:
                # panther subfamily entries have no ['signature']['entry']
                # e.g., PTHR21139:SF24
                continue

            # parse an entry, en entry is a hit vs a specific entry in a database
            motifs.append({'ipr_accession': entry['accession'],
                           'src_accession': source_accession,
                           'description':   entry['description'] or ''})

            # name = entry['name']
            # type = entry['type']

            if 'goXRefs' in entry:
                gostr = ''
                for go in entry['goXRefs']:
                    gostr += '{} ({}:{})'.format(go['id'], go['category'], go['name'])
                    if go['id'] in go_all:
                        go_all[go['id']]['source'].append(source_accession)
                    else:
                        go_all[go['id']] = {'name':   go['name'], 'category': go['category'],
                                            'source': [source_accession]}

            if 'pathwayXRefs' in entry:
                pathstr = ''
                for path in entry['pathwayXRefs']:
                    pathstr += '{} ({})'.format(path['id'], path['name'])
                    if path['databaseName'] == 'Reactome':
                        field = path['id'].split('-')
                        id = 'Reactome:{}'.format(field[2])
                    elif path['databaseName'] == 'MetaCyc':
                        id = 'Metacyc:{}'.format(path['id'])
                    elif path['databaseName'] == 'KEGG':
                        id = 'KEGG:{}'.format(path['id'])
                    else:
                        print('unknown pathway {} | {} | {}'.format(path['databaseName'],
                                                                    path['id'], path['name']))

                    if id in path_all:
                        if source_accession not in path_all[id]['source']:
                            path_all[id]['source'].append(source_accession)
                    else:
                        path_all[id] = {'name':   path['name'],
                                        'source': [source_accession]}

        return {'motifs': motifs, 'go': go_all, 'pathway': path_all}

    def run(self, show_query=False):
        """-----------------------------------------------------------------------------------------
        Construct a REST command and dispatch the job to the server
        Any previously existing jobID is overwritten

        :param show_query: boolean, print query if true
        :return: logical, True = success, False = failure
        -----------------------------------------------------------------------------------------"""
        is_success = False

        # general fields for all queries
        param = {u'email':  self.email, u'title': self.title, u'sequence': self.sequence,
                 u'output': self.output}

        if self.applications:
            # add selected applications
            param['appl'] = ','.join(self.applications)

        if self.parameters:
            for para in self.parameters:
                param[para] = self.parameters[para]

        command = self.url + 'run'
        self.response = requests.post(command, files=param, headers={'User-Agent': 'ips-client'})

        if show_query:
            # print out query if requested
            print(self.response.request.headers, '\n')
            print(self.response.request.body, '\n')

        if self.response_is_error('submitting job'):
            self.jobstatus = 'SUBMIT_ERROR'
        else:
            # success
            self.jobid = self.response.text
            self.jobstatus = 'SUBMIT_OK'
            if self.log:
                self.log_message('submitted', 'job_id={}'.format(self.jobid))

            is_success = True

        return is_success

    def status(self, log=True):
        """-----------------------------------------------------------------------------------------
        Poll job status at the server. The job is polled only once so it you want to poll
        multiple times call this method in a loop

        :return: string, status of job at server
        -----------------------------------------------------------------------------------------"""
        command = self.url + 'status/' + self.jobid
        self.response = requests.get(command)
        response_text = self.response.text.rstrip()
        if self.log > 1:
            self.log_message('polling', 'job_id={};response={}'.format(self.jobid, response_text))

        if 'FINISHED' in self.response.text:
            self.jobstatus = 'FINISHED'
            if self.log > 0:
                self.log_message('finished', 'job_id={}'.format(self.jobid))

        else:
            self.jobstatus = self.response.text
            if self.log > 0:
                self.log_message(self.jobstatus, 'job_id={}'.format(self.jobid))

        return self.jobstatus

    def result(self):
        """-----------------------------------------------------------------------------------------
        Retrieve the result

        :return: Logical True=success, False=failure
        -----------------------------------------------------------------------------------------"""
        # get the final result
        command = self.url + 'result/' + self.jobid + '/' + self.output
        self.response = requests.get(command)
        if not self.response_is_error('retrieving result'):
            # success
            self.content = self.response.text
            if self.log > 1:
                self.log_message(
                    'retrieved', 'job_id={};output_len={}'.format(self.jobid, len(self.output)))
            return True

        return False

    def response_is_error(self, task):
        """-----------------------------------------------------------------------------------------
        Return true if the response code is other than 200. Write error message to stderr if
        loglevel > 1. Task is a string describing the task that failed for inclusion in the error
        message.  The most recent response is stored in self.response

        :param task: string, text description of response being tested for error message
        :return: logical True = error, False = no error
        -----------------------------------------------------------------------------------------"""
        if self.response.status_code == 200:
            # success
            is_error = False

        else:
            # error
            is_error = True
            if self.log > 0:
                self.log_message(task, 'job_id={};status={}'.format(self.jobid,
                                                              self.response.status_code))

        return is_error

    def set_log_fh(self, fh):
        """-----------------------------------------------------------------------------------------
        The output for the log is STDERR by default.  This function allows you to change it
        -----------------------------------------------------------------------------------------"""
        self.log_fh = fh
        return fh

    def log_message(self, type, message):
        """-----------------------------------------------------------------------------------------
        write a message to the log

        types:
            not_available (for submission options)
            submitted
            polling
            finished
            retrieved
            server_error

        :type: string, type of message
        :param message: string, text of message
        :return:
        -----------------------------------------------------------------------------------------"""
        event_time = time.strftime('%d/%b/%G:%H:%M:%S', time.localtime(time.time()))
        self.log_fh.write('{}\t{}\t{}\n'.format(event_time, type, message))

        return True

    @classmethod
    def logtime(cls):
        """-----------------------------------------------------------------------------------------
        Return current time as a string. Format is 10/Oct/2000:13:55:36 which is similar to the
        common log format (without the time zone)

        DEPRECATED: The same time is generated in log_message()

        :return: string
        -----------------------------------------------------------------------------------------"""
        return time.strftime('%d/%b/%G:%H:%M:%S', time.localtime(time.time()))


def json_test():
    """---------------------------------------------------------------------------------------------
    return some json for testing

    :return:
    ---------------------------------------------------------------------------------------------"""
    return '''{
 "interproscan-version": "5.50-84.0",
"results": [ {
  "sequence" : "MAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLDPKIAVAAQNCYKVTNGAFTGEISPGMIKDCGATWVVLGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKVIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSNVSDAVAQSTRIIYGGSVTGATCKELASQPDVDGFLVGGASLKPEFVDIINAKQ",
  "md5" : "a8d44fc2c980a7677a3b54788d0fa323",
  "matches" : [ {
    "signature" : {
      "accession" : "PTHR21139",
      "name" : "TRIOSEPHOSPHATE ISOMERASE",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "PANTHER",
        "version" : "15.0"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 1,
      "end" : 249,
      "hmmStart" : 1,
      "hmmEnd" : 249,
      "hmmLength" : 249,
      "hmmBounds" : "COMPLETE",
      "envelopeStart" : 1,
      "envelopeEnd" : 249,
      "location-fragments" : [ {
        "start" : 1,
        "end" : 249,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 0.0,
    "familyName" : "Not available",
    "score" : 585.3,
    "model-ac" : "PTHR21139"
  }, {
    "signature" : {
      "accession" : "PS51440",
      "name" : "TIM_2",
      "description" : "Triosephosphate isomerase (TIM) family profile.",
      "signatureLibraryRelease" : {
        "library" : "PROSITE_PROFILES",
        "version" : "2019_11"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 6,
      "end" : 247,
      "score" : 81.047,
      "alignment" : "KFFVGGNWKMNGRKQSLGELIGTLNAaKVPAD--TEVVCAPPTAYIDFARQKLdPKIAVAAQNCYKVTNGAFTGEISPGMIKDCGATWVVLGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKVIADNVKD-WSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSNVSDAVAQSTRIIYGGSVTGATCKELASQPDVDGFLVGGASLK-PEFVDIINA",
      "location-fragments" : [ {
        "start" : 6,
        "end" : 247,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "PS51440"
  }, {
    "signature" : {
      "accession" : "MF_00147_B",
      "name" : "TIM_B",
      "description" : "Triosephosphate isomerase [tpiA].",
      "signatureLibraryRelease" : {
        "library" : "HAMAP",
        "version" : "2020_05"
      },
      "entry" : {
        "accession" : "IPR022896",
        "name" : "TrioseP_Isoase_bac/euk",
        "description" : "Triosephosphate isomerase, bacterial/eukaryotic",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        }, {
          "name" : "glycolytic process",
          "databaseName" : "GO",
          "category" : "BIOLOGICAL_PROCESS",
          "id" : "GO:0006096"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        } ]
      }
    },
    "locations" : [ {
      "start" : 5,
      "end" : 247,
      "score" : 58.311584,
      "alignment" : "Not available",
      "location-fragments" : [ {
        "start" : 5,
        "end" : 247,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "MF_00147_B"
  }, {
    "signature" : {
      "accession" : "PS00171",
      "name" : "TIM_1",
      "description" : "Triosephosphate isomerase active site.",
      "signatureLibraryRelease" : {
        "library" : "PROSITE_PATTERNS",
        "version" : "2019_11"
      },
      "entry" : {
        "accession" : "IPR020861",
        "name" : "Triosephosphate_isomerase_AS",
        "description" : "Triosephosphate isomerase, active site",
        "type" : "ACTIVE_SITE",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 164,
      "end" : 174,
      "level" : "STRONG",
      "cigarAlignment" : "11M",
      "alignment" : "AYEPVWAIGTG",
      "location-fragments" : [ {
        "start" : 164,
        "end" : 174,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "PS00171"
  }, {
    "signature" : {
      "accession" : "G3DSA:3.20.20.70",
      "name" : "Aldolase class I",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "GENE3D",
        "version" : "4.2.0"
      },
      "entry" : {
        "accession" : "IPR013785",
        "name" : "Aldolase_TIM",
        "description" : "Aldolase-type TIM barrel",
        "type" : "HOMOLOGOUS_SUPERFAMILY",
        "goXRefs" : [ {
          "name" : "catalytic activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0003824"
        } ],
        "pathwayXRefs" : [ {
          "name" : "2-chloroacrylate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7428"
        }, {
          "name" : "Geraniol degradation",
          "databaseName" : "KEGG",
          "id" : "00281"
        }, {
          "name" : "sterigmatocystin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5956"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-DME-389661"
        }, {
          "name" : "Potential therapeutics for SARS",
          "databaseName" : "Reactome",
          "id" : "R-HSA-9679191"
        }, {
          "name" : "spheroidene and spheroidenone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6286"
        }, {
          "name" : "2-deoxy-alpha-D-ribose 1-phosphate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7180"
        }, {
          "name" : "fatty acid biosynthesis initiation (type I)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5966"
        }, {
          "name" : "Platelet degranulation",
          "databaseName" : "Reactome",
          "id" : "R-MMU-114608"
        }, {
          "name" : "bacilysin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7626"
        }, {
          "name" : "Nicotinamide salvaging",
          "databaseName" : "Reactome",
          "id" : "R-DRE-197264"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-PFA-389661"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "Naphthalene degradation",
          "databaseName" : "KEGG",
          "id" : "00626"
        }, {
          "name" : "1,3,5-trimethoxybenzene biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5765"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-500753"
        }, {
          "name" : "cholesterol degradation to androstenedione I (cholesterol oxidase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6945"
        }, {
          "name" : "stellariose and mediose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6525"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-BTA-6798695"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-DME-71336"
        }, {
          "name" : "coumestrol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6332"
        }, {
          "name" : "Peroxisomal protein import",
          "databaseName" : "Reactome",
          "id" : "R-HSA-9033241"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-DRE-389661"
        }, {
          "name" : "D-galactarate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6497"
        }, {
          "name" : "3-methylbutanol biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6871"
        }, {
          "name" : "Heme biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-189451"
        }, {
          "name" : "bacimethrin and bacimethrin pyrophosphate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7564"
        }, {
          "name" : "UDP-beta-L-rhamnose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-3261"
        }, {
          "name" : "all-trans-farnesol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6859"
        }, {
          "name" : "stachyose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5337"
        }, {
          "name" : "CDP-ascarylose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5830"
        }, {
          "name" : "Histidine metabolism",
          "databaseName" : "KEGG",
          "id" : "00340"
        }, {
          "name" : "FR-900098 and FR-33289 antibiotics biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7419"
        }, {
          "name" : "dehydrophos biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6682"
        }, {
          "name" : "superpathway of mycolate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6113"
        }, {
          "name" : "D-galacturonate degradation III",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6491"
        }, {
          "name" : "podophyllotoxin glucosides metabolism",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7491"
        }, {
          "name" : "bile acid 7alpha-dehydroxylation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7754"
        }, {
          "name" : "raspberry ketone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5393"
        }, {
          "name" : "Sialic acid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-DRE-4085001"
        }, {
          "name" : "7-(3-amino-3-carboxypropyl)-wyosine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7286"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "type I lipoteichoic acid biosynthesis (S. aureus)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7817"
        }, {
          "name" : "Glyoxylate and dicarboxylate metabolism",
          "databaseName" : "KEGG",
          "id" : "00630"
        }, {
          "name" : "vindoline, vindorosine and vinblastine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5292"
        }, {
          "name" : "thiamine diphosphate biosynthesis I (E. coli)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6894"
        }, {
          "name" : "glyphosate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7806"
        }, {
          "name" : "thymine degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6430"
        }, {
          "name" : "methylerythritol phosphate pathway II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7560"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-RNO-6798695"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-XTR-389661"
        }, {
          "name" : "jasmonic acid biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-735"
        }, {
          "name" : "CMP-N-acetylneuraminate biosynthesis I (eukaryotes)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6138"
        }, {
          "name" : "detoxification of reactive carbonyls in chloroplasts",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6786"
        }, {
          "name" : "Peroxisomal lipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-DDI-390918"
        }, {
          "name" : "Sialic acid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-MMU-4085001"
        }, {
          "name" : "equisetin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7669"
        }, {
          "name" : "Lipopolysaccharide biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00540"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-SPO-6798695"
        }, {
          "name" : "Nitrogen metabolism",
          "databaseName" : "KEGG",
          "id" : "00910"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-SCE-6798695"
        }, {
          "name" : "bacterial bioluminescence",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7723"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "astaxanthin biosynthesis (flowering plants)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7636"
        }, {
          "name" : "CMP-8-amino-3,8-dideoxy-D-manno-octulosonate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7674"
        }, {
          "name" : "(S)-lactate fermentation to propanoate, acetate and hydrogen",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8086"
        }, {
          "name" : "3-dehydroquinate biosynthesis II (archaea)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6160"
        }, {
          "name" : "Pantothenate and CoA biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00770"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "hesperitin glycoside biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5105"
        }, {
          "name" : "bixin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5305"
        }, {
          "name" : "dTDP-beta-L-mycarose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6976"
        }, {
          "name" : "Pyrimidine catabolism",
          "databaseName" : "Reactome",
          "id" : "R-DRE-73621"
        }, {
          "name" : "quinate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6416"
        }, {
          "name" : "(+)-pisatin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2467"
        }, {
          "name" : "choline biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-3542"
        }, {
          "name" : "dTDP-beta-L-digitoxose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7657"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-PFA-6798695"
        }, {
          "name" : "Nicotinate metabolism",
          "databaseName" : "Reactome",
          "id" : "R-MMU-196807"
        }, {
          "name" : "Nicotinamide salvaging",
          "databaseName" : "Reactome",
          "id" : "R-MMU-197264"
        }, {
          "name" : "4-nitrobenzoate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2381"
        }, {
          "name" : "aclacinomycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7354"
        }, {
          "name" : "Galactose metabolism",
          "databaseName" : "KEGG",
          "id" : "00052"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70350"
        }, {
          "name" : "tetrapyrrole biosynthesis II (from glycine)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5189"
        }, {
          "name" : "toyocamycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6720"
        }, {
          "name" : "linezolid resistance",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6828"
        }, {
          "name" : "3-hydroxy-4-methyl-anthranilate biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7717"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70350"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-HSA-6798695"
        }, {
          "name" : "blasticidin S biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7570"
        }, {
          "name" : "Purine salvage",
          "databaseName" : "Reactome",
          "id" : "R-HSA-74217"
        }, {
          "name" : "thiamine diphosphate biosynthesis II (Bacillus)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6893"
        }, {
          "name" : "salicortin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6763"
        }, {
          "name" : "lychnose and isolychnose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6524"
        }, {
          "name" : "photorespiration",
          "databaseName" : "MetaCyc",
          "id" : "PWY-181"
        }, {
          "name" : "aerobic respiration II (cytochrome c) (yeast)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7279"
        }, {
          "name" : "spermidine hydroxycinnamic acid conjugates biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6442"
        }, {
          "name" : "homogalacturonan biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1061"
        }, {
          "name" : "UDP-N-acetyl-beta-L-quinovosamine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7331"
        }, {
          "name" : "4-amino-2-methyl-5-diphosphomethylpyrimidine biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6890"
        }, {
          "name" : "phenylethanol glycoconjugate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7074"
        }, {
          "name" : "coniferyl alcohol 9-methyl ester biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7643"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-MMU-389661"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Nicotinamide salvaging",
          "databaseName" : "Reactome",
          "id" : "R-HSA-197264"
        }, {
          "name" : "Peroxisomal lipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-RNO-390918"
        }, {
          "name" : "CMP-legionaminate biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6749"
        }, {
          "name" : "4,4'-disulfanediyldibutanoate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7898"
        }, {
          "name" : "oleandomycin activation/inactivation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6972"
        }, {
          "name" : "Nicotinate metabolism",
          "databaseName" : "Reactome",
          "id" : "R-DDI-196807"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-73817"
        }, {
          "name" : "phenylpropanoids methylation (ice plant)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7498"
        }, {
          "name" : "thiamine diphosphate salvage II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6897"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "benzoyl-CoA biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6458"
        }, {
          "name" : "cyclohexane-1-carboxyl-CoA biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8044"
        }, {
          "name" : "Glucosinolate biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00966"
        }, {
          "name" : "kavain biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8103"
        }, {
          "name" : "justicidin B biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6824"
        }, {
          "name" : "ephedrine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5883"
        }, {
          "name" : "Glycosphingolipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-DDI-1660662"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "streptovaricin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8035"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-SPO-71336"
        }, {
          "name" : "autoinducer AI-2 biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6153"
        }, {
          "name" : "2-deoxy-D-ribose degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8060"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "Pyrimidine metabolism",
          "databaseName" : "KEGG",
          "id" : "00240"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70350"
        }, {
          "name" : "anditomin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7599"
        }, {
          "name" : "xylan biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5800"
        }, {
          "name" : "Glycosphingolipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-MMU-1660662"
        }, {
          "name" : "jadomycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6679"
        }, {
          "name" : "Pentose and glucuronate interconversions",
          "databaseName" : "KEGG",
          "id" : "00040"
        }, {
          "name" : "Insulin effects increased synthesis of Xylulose-5-Phosphate",
          "databaseName" : "Reactome",
          "id" : "R-MMU-163754"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-DDI-6798695"
        }, {
          "name" : "L-glucose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7130"
        }, {
          "name" : "lincomycin A biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6955"
        }, {
          "name" : "preQ0 biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6703"
        }, {
          "name" : "acridone alkaloid biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5958"
        }, {
          "name" : "Biosynthesis of vancomycin group antibiotics",
          "databaseName" : "KEGG",
          "id" : "01055"
        }, {
          "name" : "sitosterol degradation to androstenedione",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6948"
        }, {
          "name" : "adenosine nucleotides degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6596"
        }, {
          "name" : "Nicotinate and nicotinamide metabolism",
          "databaseName" : "KEGG",
          "id" : "00760"
        }, {
          "name" : "N-acetyl-D-galactosamine degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7077"
        }, {
          "name" : "Biosynthesis of various secondary metabolites - part 3",
          "databaseName" : "KEGG",
          "id" : "00997"
        }, {
          "name" : "3-hydroxy-L-homotyrosine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7649"
        }, {
          "name" : "1,4-dihydroxy-6-naphthoate biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7374"
        }, {
          "name" : "6-gingerol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6835"
        }, {
          "name" : "L-lysine biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2941"
        }, {
          "name" : "Amino sugar and nucleotide sugar metabolism",
          "databaseName" : "KEGG",
          "id" : "00520"
        }, {
          "name" : "dTDP-alpha-D-olivose, dTDP-alpha-D-oliose and dTDP-alpha-D-mycarose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6973"
        }, {
          "name" : "Methanobacterium thermoautotrophicum biosynthetic metabolism",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6146"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "heme b biosynthesis III (from siroheme)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7552"
        }, {
          "name" : "dTDP-beta-L-4-epi-vancosamine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7440"
        }, {
          "name" : "shinorine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7751"
        }, {
          "name" : "citronellol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6670"
        }, {
          "name" : "stearate biosynthesis I (animals)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5972"
        }, {
          "name" : "sorgoleone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5987"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "3,6-anhydro-alpha-L-galactopyranose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7562"
        }, {
          "name" : "bisabolene biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7102"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-500753"
        }, {
          "name" : "Sialic acid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-HSA-4085001"
        }, {
          "name" : "anteiso-branched-chain fatty acid biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8173"
        }, {
          "name" : "yangonin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8104"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "ulvan degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7647"
        }, {
          "name" : "Lysine degradation",
          "databaseName" : "KEGG",
          "id" : "00310"
        }, {
          "name" : "Peroxisomal lipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-BTA-390918"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-CEL-71336"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "odd iso-branched-chain fatty acid biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8174"
        }, {
          "name" : "Nicotinamide salvaging",
          "databaseName" : "Reactome",
          "id" : "R-RNO-197264"
        }, {
          "name" : "mupirocin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8012"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-DME-6798695"
        }, {
          "name" : "factor 420 biosynthesis I (archaea)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8112"
        }, {
          "name" : "papaverine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7363"
        }, {
          "name" : "daphnin interconversion",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7056"
        }, {
          "name" : "Biotin transport and metabolism",
          "databaseName" : "Reactome",
          "id" : "R-CEL-196780"
        }, {
          "name" : "Porphyrin and chlorophyll metabolism",
          "databaseName" : "KEGG",
          "id" : "00860"
        }, {
          "name" : "Synthesis of Ketone Bodies",
          "databaseName" : "Reactome",
          "id" : "R-DRE-77111"
        }, {
          "name" : "Biotin transport and metabolism",
          "databaseName" : "Reactome",
          "id" : "R-SPO-196780"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-73817"
        }, {
          "name" : "CS/DS degradation",
          "databaseName" : "Reactome",
          "id" : "R-MMU-2024101"
        }, {
          "name" : "Pyrimidine catabolism",
          "databaseName" : "Reactome",
          "id" : "R-RNO-73621"
        }, {
          "name" : "dalcochinin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5821"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "pederin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8049"
        }, {
          "name" : "Valine, leucine and isoleucine degradation",
          "databaseName" : "KEGG",
          "id" : "00280"
        }, {
          "name" : "grixazone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7153"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-SCE-389661"
        }, {
          "name" : "Ubiquinone and other terpenoid-quinone biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00130"
        }, {
          "name" : "sangivamycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6721"
        }, {
          "name" : "Arginine biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00220"
        }, {
          "name" : "Nicotinate metabolism",
          "databaseName" : "Reactome",
          "id" : "R-SCE-196807"
        }, {
          "name" : "bile acid 7beta-dehydroxylation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8134"
        }, {
          "name" : "Entner-Doudoroff pathway III (semi-phosphorylative)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2221"
        }, {
          "name" : "thiamine phosphate formation from pyrithiamine and oxythiamine (yeast)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7357"
        }, {
          "name" : "echinomycin and triostin A biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7735"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "gamma-coniciene and coniine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5748"
        }, {
          "name" : "mevalonate pathway II (haloarchaea)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6174"
        }, {
          "name" : "gallate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6707"
        }, {
          "name" : "Glycine, serine and threonine metabolism",
          "databaseName" : "KEGG",
          "id" : "00260"
        }, {
          "name" : "Glycosphingolipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-HSA-1660662"
        }, {
          "name" : "methanofuran biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5254"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "KEGG",
          "id" : "00030"
        }, {
          "name" : "8-O-methylfusarubin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7689"
        }, {
          "name" : "benzoate biosynthesis III (CoA-dependent, non-beta-oxidative)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6446"
        }, {
          "name" : "mithramycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7045"
        }, {
          "name" : "Sialic acid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-XTR-4085001"
        }, {
          "name" : "Lysine biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00300"
        }, {
          "name" : "dTDP-L-daunosamine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7814"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "paromamine biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7022"
        }, {
          "name" : "L-isoleucine biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5101"
        }, {
          "name" : "Nicotinate metabolism",
          "databaseName" : "Reactome",
          "id" : "R-RNO-196807"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "dipicolinate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8088"
        }, {
          "name" : "colchicine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5846"
        }, {
          "name" : "gentamicin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7025"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70350"
        }, {
          "name" : "dTDP-beta-L-evernitrose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7439"
        }, {
          "name" : "benzoate biosynthesis I (CoA-dependent, beta-oxidative)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6443"
        }, {
          "name" : "ecdysone and 20-hydroxyecdysone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7300"
        }, {
          "name" : "mycofactocin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8122"
        }, {
          "name" : "aurofusarin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7695"
        }, {
          "name" : "gadusol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7752"
        }, {
          "name" : "Arginine and proline metabolism",
          "databaseName" : "KEGG",
          "id" : "00330"
        }, {
          "name" : "vestitol and sativan biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5729"
        }, {
          "name" : "pinitol biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6739"
        }, {
          "name" : "CS/DS degradation",
          "databaseName" : "Reactome",
          "id" : "R-RNO-2024101"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70350"
        }, {
          "name" : "CS/DS degradation",
          "databaseName" : "Reactome",
          "id" : "R-HSA-2024101"
        }, {
          "name" : "3-amino-5-hydroxybenzoate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5979"
        }, {
          "name" : "apicidin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7667"
        }, {
          "name" : "aucuparin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7739"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-500753"
        }, {
          "name" : "Pyrimidine catabolism",
          "databaseName" : "Reactome",
          "id" : "R-MMU-73621"
        }, {
          "name" : "5-methoxy-6-methylbenzimidazole biosynthesis (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8095"
        }, {
          "name" : "L-arabinose degradation III",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5517"
        }, {
          "name" : "Nicotinate metabolism",
          "databaseName" : "Reactome",
          "id" : "R-HSA-196807"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "okenone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7591"
        }, {
          "name" : "tRNA modification in the nucleus and cytosol",
          "databaseName" : "Reactome",
          "id" : "R-HSA-6782315"
        }, {
          "name" : "propanethial S-oxide biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5707"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-73817"
        }, {
          "name" : "Insulin effects increased synthesis of Xylulose-5-Phosphate",
          "databaseName" : "Reactome",
          "id" : "R-DME-163754"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Peroxisomal protein import",
          "databaseName" : "Reactome",
          "id" : "R-MMU-9033241"
        }, {
          "name" : "pyrroloquinoline quinone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6420"
        }, {
          "name" : "4-nitrotoluene degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5644"
        }, {
          "name" : "TALDO1 deficiency: failed conversion of SH7P, GA3P to Fru(6)P, E4P",
          "databaseName" : "Reactome",
          "id" : "R-HSA-6791055"
        }, {
          "name" : "UMP biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5686"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "dechlorogriseofulvin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7655"
        }, {
          "name" : "CMP-N-acetyl-7-O-acetylneuraminate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7529"
        }, {
          "name" : "Phosphonate and phosphinate metabolism",
          "databaseName" : "KEGG",
          "id" : "00440"
        }, {
          "name" : "NAD salvage (plants)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5381"
        }, {
          "name" : "ginkgotoxin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8077"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-500753"
        }, {
          "name" : "bikaverin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7692"
        }, {
          "name" : "L-tyrosine biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-3461"
        }, {
          "name" : "Synthesis of Ketone Bodies",
          "databaseName" : "Reactome",
          "id" : "R-RNO-77111"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        }, {
          "name" : "syringate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6339"
        }, {
          "name" : "even iso-branched-chain fatty acid biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8175"
        }, {
          "name" : "mono-trans, poly-cis decaprenyl phosphate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6383"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Nicotinate metabolism",
          "databaseName" : "Reactome",
          "id" : "R-SSC-196807"
        }, {
          "name" : "methyl-coenzyme M oxidation to CO2",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5209"
        }, {
          "name" : "mevalonate pathway I (eukaryotes and bacteria)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-922"
        }, {
          "name" : "Biotin transport and metabolism",
          "databaseName" : "Reactome",
          "id" : "R-MMU-196780"
        }, {
          "name" : "pentose phosphate pathway (oxidative branch) II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7796"
        }, {
          "name" : "sulfur volatiles biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6736"
        }, {
          "name" : "molybdopterin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6823"
        }, {
          "name" : "pinobanksin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5059"
        }, {
          "name" : "geodin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7079"
        }, {
          "name" : "Phenazine biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00405"
        }, {
          "name" : "Peroxisomal protein import",
          "databaseName" : "Reactome",
          "id" : "R-DDI-9033241"
        }, {
          "name" : "ergosterol biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7154"
        }, {
          "name" : "Sphingolipid metabolism",
          "databaseName" : "KEGG",
          "id" : "00600"
        }, {
          "name" : "daphnetin modification",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7055"
        }, {
          "name" : "divinyl ether biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5409"
        }, {
          "name" : "Rubisco shunt",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5723"
        }, {
          "name" : "Platelet degranulation",
          "databaseName" : "Reactome",
          "id" : "R-PFA-114608"
        }, {
          "name" : "Interferon alpha/beta signaling",
          "databaseName" : "Reactome",
          "id" : "R-HSA-909733"
        }, {
          "name" : "hopanoid biosynthesis (bacteria)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7072"
        }, {
          "name" : "Heme biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-189451"
        }, {
          "name" : "Benzoate degradation",
          "databaseName" : "KEGG",
          "id" : "00362"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-DRE-6798695"
        }, {
          "name" : "diacylglyceryl-N,N,N-trimethylhomoserine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6795"
        }, {
          "name" : "digitoxigenin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6032"
        }, {
          "name" : "Lipoic acid metabolism",
          "databaseName" : "KEGG",
          "id" : "00785"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Citrate cycle (TCA cycle)",
          "databaseName" : "KEGG",
          "id" : "00020"
        }, {
          "name" : "polymethylated quercetin glucoside biosynthesis II - quercetagetin series (Chrysosplenium)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7151"
        }, {
          "name" : "rhizocticin A and B biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7510"
        }, {
          "name" : "4-hydroxy-2(1H)-quinolone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6661"
        }, {
          "name" : "chaxamycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8036"
        }, {
          "name" : "D-altritol and galactitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7862"
        }, {
          "name" : "fatty acid biosynthesis initiation (type II)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-4381"
        }, {
          "name" : "traumatin and (Z)-3-hexen-1-yl acetate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5410"
        }, {
          "name" : "Purine salvage",
          "databaseName" : "Reactome",
          "id" : "R-CEL-74217"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Heme biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-189451"
        }, {
          "name" : "griseofulvin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7653"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "dTDP-4-O-demethyl-beta-L-noviose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7301"
        }, {
          "name" : "dTDP-beta-L-megosamine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7104"
        }, {
          "name" : "UMP biosynthesis III",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7791"
        }, {
          "name" : "daidzin and daidzein degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6996"
        }, {
          "name" : "indolmycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7770"
        }, {
          "name" : "triclosan resistance",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7096"
        }, {
          "name" : "Folate biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00790"
        }, {
          "name" : "L-glutamine biosynthesis III",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6549"
        }, {
          "name" : "Sialic acid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-DME-4085001"
        }, {
          "name" : "brassicicene C biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7517"
        }, {
          "name" : "proline betaine degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8109"
        }, {
          "name" : "neomycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7016"
        }, {
          "name" : "viridicatumtoxin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7659"
        }, {
          "name" : "3-chlorobenzoate degradation I (via chlorocatechol)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6088"
        }, {
          "name" : "bryostatin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8047"
        }, {
          "name" : "bisphenol A degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7757"
        }, {
          "name" : "Interaction With Cumulus Cells And The Zona Pellucida",
          "databaseName" : "Reactome",
          "id" : "R-HSA-2534343"
        }, {
          "name" : "4-deoxy-L-threo-hex-4-enopyranuronate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6507"
        }, {
          "name" : "D-glucuronate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6501"
        }, {
          "name" : "thiazole component of thiamine diphosphate biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6891"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Glycerophospholipid metabolism",
          "databaseName" : "KEGG",
          "id" : "00564"
        }, {
          "name" : "methylthiopropanonate degradation II (demethylation)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6045"
        }, {
          "name" : "3PG-factor 420 biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8113"
        }, {
          "name" : "daidzein conjugates interconversion",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2343"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-73817"
        }, {
          "name" : "Platelet degranulation",
          "databaseName" : "Reactome",
          "id" : "R-DDI-114608"
        }, {
          "name" : "Heme biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-189451"
        }, {
          "name" : "CDP-6-deoxy-D-gulose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8139"
        }, {
          "name" : "rutin degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6848"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Phenylalanine metabolism",
          "databaseName" : "KEGG",
          "id" : "00360"
        }, {
          "name" : "superpathway of pterocarpan biosynthesis (via formononetin)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2229"
        }, {
          "name" : "dimycocerosyl triglycosyl phenolphthiocerol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7743"
        }, {
          "name" : "Alanine, aspartate and glutamate metabolism",
          "databaseName" : "KEGG",
          "id" : "00250"
        }, {
          "name" : "FeMo cofactor biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7710"
        }, {
          "name" : "Gene and protein expression by JAK-STAT signaling after Interleukin-12 stimulation",
          "databaseName" : "Reactome",
          "id" : "R-HSA-8950505"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-500753"
        }, {
          "name" : "fusarin C biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7673"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-XTR-6798695"
        }, {
          "name" : "MPS IX - Natowicz syndrome",
          "databaseName" : "Reactome",
          "id" : "R-HSA-2206280"
        }, {
          "name" : "platensimycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8179"
        }, {
          "name" : "Xylene degradation",
          "databaseName" : "KEGG",
          "id" : "00622"
        }, {
          "name" : "arctigenin and isoarctigenin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7631"
        }, {
          "name" : "Benzoxazinoid biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00402"
        }, {
          "name" : "Glycerolipid metabolism",
          "databaseName" : "KEGG",
          "id" : "00561"
        }, {
          "name" : "indole-3-acetate activation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1921"
        }, {
          "name" : "N-acetylneuraminate and N-acetylmannosamine degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7581"
        }, {
          "name" : "alpha-Linolenic acid metabolism",
          "databaseName" : "KEGG",
          "id" : "00592"
        }, {
          "name" : "lipoate biosynthesis and incorporation III (Bacillus)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6987"
        }, {
          "name" : "cis-geranyl-CoA degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6672"
        }, {
          "name" : "emetine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7135"
        }, {
          "name" : "mycobacterial sulfolipid biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7746"
        }, {
          "name" : "salinosporamide A biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6627"
        }, {
          "name" : "protein N-glycosylation (Haloferax volcanii)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7661"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-73817"
        }, {
          "name" : "luteolin triglucuronide degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7445"
        }, {
          "name" : "Thiamine metabolism",
          "databaseName" : "KEGG",
          "id" : "00730"
        }, {
          "name" : "L-lysine biosynthesis VI",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5097"
        }, {
          "name" : "Glycosphingolipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-XTR-1660662"
        }, {
          "name" : "hypericin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5780"
        }, {
          "name" : "oleate biosynthesis I (plants)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5147"
        }, {
          "name" : "pyridoxal 5'-phosphate biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6466"
        }, {
          "name" : "Defective HLCS causes multiple carboxylase deficiency",
          "databaseName" : "Reactome",
          "id" : "R-HSA-3371599"
        }, {
          "name" : "spirilloxanthin and 2,2'-diketo-spirilloxanthin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6581"
        }, {
          "name" : "thiocoraline biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7737"
        }, {
          "name" : "seleno-amino acid detoxification and volatilization II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6935"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "inosine 5'-phosphate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5695"
        }, {
          "name" : "polymethylated quercetin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7161"
        }, {
          "name" : "palmitate biosynthesis (type II fatty acid synthase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5971"
        }, {
          "name" : "Purine salvage",
          "databaseName" : "Reactome",
          "id" : "R-RNO-74217"
        }, {
          "name" : "fumigaclavine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7059"
        }, {
          "name" : "brassinosteroid biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2582"
        }, {
          "name" : "eriodictyol C-glucosylation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7188"
        }, {
          "name" : "3-methylquinoline degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-721"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-SCE-71336"
        }, {
          "name" : "Pyruvate metabolism",
          "databaseName" : "KEGG",
          "id" : "00620"
        }, {
          "name" : "D-galacturonate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6486"
        }, {
          "name" : "Purine salvage",
          "databaseName" : "Reactome",
          "id" : "R-MMU-74217"
        }, {
          "name" : "Dioxin degradation",
          "databaseName" : "KEGG",
          "id" : "00621"
        }, {
          "name" : "Nicotinamide salvaging",
          "databaseName" : "Reactome",
          "id" : "R-CEL-197264"
        }, {
          "name" : "Fatty acid biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00061"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "thiamine diphosphate salvage IV (yeast)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7356"
        }, {
          "name" : "DIBOA-glucoside biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6949"
        }, {
          "name" : "thiazole component of thiamine diphosphate biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6892"
        }, {
          "name" : "Biotin transport and metabolism",
          "databaseName" : "Reactome",
          "id" : "R-HSA-196780"
        }, {
          "name" : "Insulin effects increased synthesis of Xylulose-5-Phosphate",
          "databaseName" : "Reactome",
          "id" : "R-DDI-163754"
        }, {
          "name" : "tunicamycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7821"
        }, {
          "name" : "vanillin biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5665"
        }, {
          "name" : "rutin degradation (plants)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7134"
        }, {
          "name" : "Biotin transport and metabolism",
          "databaseName" : "Reactome",
          "id" : "R-RNO-196780"
        }, {
          "name" : "Sialic acid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-SSC-4085001"
        }, {
          "name" : "Platelet degranulation",
          "databaseName" : "Reactome",
          "id" : "R-HSA-114608"
        }, {
          "name" : "chitin degradation I (archaea)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6855"
        }, {
          "name" : "oleate biosynthesis IV (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7664"
        }, {
          "name" : "androsrtendione degradation II (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8152"
        }, {
          "name" : "9-lipoxygenase and 9-hydroperoxide lyase pathway",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5408"
        }, {
          "name" : "Monobactam biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00261"
        }, {
          "name" : "L-phenylalanine biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-3462"
        }, {
          "name" : "mevalonate pathway III (Thermoplasma)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7524"
        }, {
          "name" : "brassinosteroid biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-699"
        }, {
          "name" : "beta-Alanine metabolism",
          "databaseName" : "KEGG",
          "id" : "00410"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "L-dopa degradation II (bacterial)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8110"
        }, {
          "name" : "alkane biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7033"
        }, {
          "name" : "dTDP-beta-L-olivose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6974"
        }, {
          "name" : "3,5-dimethoxytoluene biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7076"
        }, {
          "name" : "7-dehydroporiferasterol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7155"
        }, {
          "name" : "stachyose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6527"
        }, {
          "name" : "Pyrimidine catabolism",
          "databaseName" : "Reactome",
          "id" : "R-DDI-73621"
        }, {
          "name" : "phthiocerol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7741"
        }, {
          "name" : "artemisinin and arteannuin B biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5195"
        }, {
          "name" : "polymethylated quercetin glucoside biosynthesis I - quercetin series (Chrysosplenium)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7150"
        }, {
          "name" : "CMP-legionaminate biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7131"
        }, {
          "name" : "alpha-dystroglycan glycosylation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7981"
        }, {
          "name" : "Sialic acid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-RNO-4085001"
        }, {
          "name" : "GDP-6-deoxy-D-altro-heptose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7610"
        }, {
          "name" : "1,4-dihydroxy-6-naphthoate biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7371"
        }, {
          "name" : "androstenedione degradation I (aerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6944"
        }, {
          "name" : "methylglyoxal degradation V",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5458"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70350"
        }, {
          "name" : "guanosine ribonucleotides de novo biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7221"
        }, {
          "name" : "palmitoleate biosynthesis II (plants and bacteria)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5366"
        }, {
          "name" : "(5Z)-dodecenoate biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7858"
        }, {
          "name" : "superpathway of benzoxazinoid glucosides biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-4161"
        }, {
          "name" : "cholesterol degradation to androstenedione II (cholesterol dehydrogenase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6946"
        }, {
          "name" : "phytosterol biosynthesis (plants)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2541"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-CEL-6798695"
        }, {
          "name" : "chitin derivatives degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6906"
        }, {
          "name" : "L-ascorbate degradation II (bacterial, aerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6961"
        }, {
          "name" : "staurosporine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6346"
        }, {
          "name" : "GDP-6-deoxy-D-manno-heptose biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7613"
        }, {
          "name" : "Insulin effects increased synthesis of Xylulose-5-Phosphate",
          "databaseName" : "Reactome",
          "id" : "R-HSA-163754"
        }, {
          "name" : "p-HBAD biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7745"
        }, {
          "name" : "BMAL1:CLOCK,NPAS2 activates circadian gene expression",
          "databaseName" : "Reactome",
          "id" : "R-HSA-1368108"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-500753"
        }, {
          "name" : "Secondary bile acid biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00121"
        }, {
          "name" : "Neutrophil degranulation",
          "databaseName" : "Reactome",
          "id" : "R-MMU-6798695"
        }, {
          "name" : "Nicotinamide salvaging",
          "databaseName" : "Reactome",
          "id" : "R-DME-197264"
        }, {
          "name" : "formaldehyde oxidation VI (H4MPT pathway)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1723"
        }, {
          "name" : "epiberberine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8031"
        }, {
          "name" : "Hyaluronan uptake and degradation",
          "databaseName" : "Reactome",
          "id" : "R-RNO-2160916"
        }, {
          "name" : "Chorismate via Shikimate Pathway",
          "databaseName" : "Reactome",
          "id" : "R-MTU-964903"
        }, {
          "name" : "gliotoxin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7533"
        }, {
          "name" : "plastoquinol-9 biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6978"
        }, {
          "name" : "Interaction With Cumulus Cells And The Zona Pellucida",
          "databaseName" : "Reactome",
          "id" : "R-MMU-2534343"
        }, {
          "name" : "tetrapyrrole biosynthesis I (from glutamate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5188"
        }, {
          "name" : "L-homomethionine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1186"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-HSA-71336"
        }, {
          "name" : "UDP-N-acetyl-alpha-D-galactosaminuronate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7336"
        }, {
          "name" : "isoprene biosynthesis II (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7391"
        }, {
          "name" : "Molybdenum cofactor biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-947581"
        }, {
          "name" : "D-Glutamine and D-glutamate metabolism",
          "databaseName" : "KEGG",
          "id" : "00471"
        }, {
          "name" : "ammonia assimilation cycle I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6963"
        }, {
          "name" : "furaneol and mesifurane biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5975"
        }, {
          "name" : "tetracenomycin C biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7485"
        }, {
          "name" : "fructan degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-862"
        }, {
          "name" : "yatein biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7629"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70350"
        }, {
          "name" : "petroselinate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5367"
        }, {
          "name" : "CMP-pseudaminate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6143"
        }, {
          "name" : "eupatolitin 3-O-glucoside biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7157"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-73817"
        }, {
          "name" : "superpathway of scopolin and esculin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7186"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-BTA-389661"
        }, {
          "name" : "octanoyl-[acyl-carrier protein] biosynthesis (mitochondria, yeast)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7388"
        }, {
          "name" : "sesamin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5469"
        }, {
          "name" : "superpathway of nicotine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7342"
        }, {
          "name" : "Terpenoid backbone biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00900"
        }, {
          "name" : "Insulin effects increased synthesis of Xylulose-5-Phosphate",
          "databaseName" : "Reactome",
          "id" : "R-SPO-163754"
        }, {
          "name" : "beta-alanine betaine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-4021"
        }, {
          "name" : "Drug metabolism - other enzymes",
          "databaseName" : "KEGG",
          "id" : "00983"
        }, {
          "name" : "esculetin modification",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7058"
        }, {
          "name" : "quebrachitol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8078"
        }, {
          "name" : "validamycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5818"
        }, {
          "name" : "mevalonate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5074"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "isoflavonoid biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2083"
        }, {
          "name" : "arginomycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7569"
        }, {
          "name" : "8-O-methylated benzoxazinoid glucoside biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8076"
        }, {
          "name" : "elloramycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7483"
        }, {
          "name" : "factor 420 biosynthesis II (mycobacteria)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5198"
        }, {
          "name" : "polymethylated kaempferol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7163"
        }, {
          "name" : "Molybdenum cofactor biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-947581"
        }, {
          "name" : "Purine metabolism",
          "databaseName" : "KEGG",
          "id" : "00230"
        }, {
          "name" : "glucosinolate biosynthesis from tryptophan",
          "databaseName" : "MetaCyc",
          "id" : "PWY-601"
        }, {
          "name" : "trans, trans-farnesyl diphosphate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5123"
        }, {
          "name" : "methylwyosine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7285"
        }, {
          "name" : "Heme biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-189451"
        }, {
          "name" : "Vitamin B6 metabolism",
          "databaseName" : "KEGG",
          "id" : "00750"
        }, {
          "name" : "Glycosphingolipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-CEL-1660662"
        }, {
          "name" : "Biotin transport and metabolism",
          "databaseName" : "Reactome",
          "id" : "R-SCE-196780"
        }, {
          "name" : "starch degradation IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6735"
        }, {
          "name" : "chorismate biosynthesis from 3-dehydroquinate",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6163"
        }, {
          "name" : "Molybdenum cofactor biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-947581"
        }, {
          "name" : "alkylnitronates degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-723"
        }, {
          "name" : "tetrahydromethanopterin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6148"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "Methane metabolism",
          "databaseName" : "KEGG",
          "id" : "00680"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-HSA-389661"
        }, {
          "name" : "kappa-carrageenan degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6821"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-CEL-389661"
        }, {
          "name" : "palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6282"
        }, {
          "name" : "Heme biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-BTA-189451"
        }, {
          "name" : "naphthalene degradation (aerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5427"
        }, {
          "name" : "dimethylsulfoniopropanoate biosynthesis II (Spartina)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6055"
        }, {
          "name" : "triethylamine degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7085"
        }, {
          "name" : "paromomycin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7018"
        }, {
          "name" : "toluene degradation to benzoyl-CoA (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-81"
        }, {
          "name" : "ammonia assimilation cycle II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6964"
        }, {
          "name" : "phosphinothricin tripeptide biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6322"
        }, {
          "name" : "6-methylpretetramide biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7811"
        }, {
          "name" : "daunorubicin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7352"
        }, {
          "name" : "flaviolin dimer and mompain biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7513"
        }, {
          "name" : "Valine, leucine and isoleucine biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00290"
        }, {
          "name" : "gibberellin inactivation II (methylation)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6477"
        }, {
          "name" : "10,13-epoxy-11-methyl-octadecadienoate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7691"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-SPO-389661"
        }, {
          "name" : "limonene degradation IV (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8029"
        }, {
          "name" : "quercetin glucoside degradation (Allium)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7133"
        }, {
          "name" : "Phenylalanine, tyrosine and tryptophan biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00400"
        }, {
          "name" : "Insulin effects increased synthesis of Xylulose-5-Phosphate",
          "databaseName" : "Reactome",
          "id" : "R-SCE-163754"
        }, {
          "name" : "Polycyclic aromatic hydrocarbon degradation",
          "databaseName" : "KEGG",
          "id" : "00624"
        }, {
          "name" : "Heme biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-189451"
        }, {
          "name" : "Platelet degranulation",
          "databaseName" : "Reactome",
          "id" : "R-CEL-114608"
        }, {
          "name" : "pentose phosphate pathway (non-oxidative branch) II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8178"
        }, {
          "name" : "pinocembrin C-glucosylation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7189"
        }, {
          "name" : "Neomycin, kanamycin and gentamicin biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00524"
        }, {
          "name" : "Synthesis of wybutosine at G37 of tRNA(Phe)",
          "databaseName" : "Reactome",
          "id" : "R-HSA-6782861"
        }, {
          "name" : "ansatrienin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8040"
        }, {
          "name" : "rot-2'-enonate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6427"
        }, {
          "name" : "ajmaline and sarpagine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5301"
        }, {
          "name" : "D-fructuronate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7242"
        }, {
          "name" : "Glycosphingolipid biosynthesis - globo and isoglobo series",
          "databaseName" : "KEGG",
          "id" : "00603"
        }, {
          "name" : "cyanidin diglucoside biosynthesis (acyl-glucose dependent)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7256"
        }, {
          "name" : "mitochondrial NADPH production (yeast)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7269"
        }, {
          "name" : "phenolphthiocerol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7742"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-MMU-71336"
        }, {
          "name" : "uracil degradation I (reductive)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-3982"
        }, {
          "name" : "thiamine diphosphate biosynthesis IV (eukaryotes)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6908"
        }, {
          "name" : "phosalacine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7769"
        }, {
          "name" : "4-amino-2-methyl-5-diphosphomethylpyrimidine biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7282"
        }, {
          "name" : "NAD biosynthesis from 2-amino-3-carboxymuconate semialdehyde",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5653"
        }, {
          "name" : "D-glucarate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6499"
        }, {
          "name" : "D-galactosamine and N-acetyl-D-galactosamine degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7395"
        }, {
          "name" : "salvigenin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7325"
        }, {
          "name" : "5-methoxybenzimidazole biosynthesis (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8096"
        }, {
          "name" : "reductive TCA cycle II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5392"
        }, {
          "name" : "Aminobenzoate degradation",
          "databaseName" : "KEGG",
          "id" : "00627"
        }, {
          "name" : "formaldehyde assimilation II (assimilatory RuMP Cycle)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1861"
        }, {
          "name" : "glutaminyl-tRNAgln biosynthesis via transamidation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5921"
        }, {
          "name" : "dTDP-D-desosamine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6942"
        }, {
          "name" : "Platelet degranulation",
          "databaseName" : "Reactome",
          "id" : "R-DME-114608"
        }, {
          "name" : "Butanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00650"
        }, {
          "name" : "butanol and isobutanol biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7396"
        }, {
          "name" : "L-tyrosine biosynthesis III",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6120"
        }, {
          "name" : "gondoate biosynthesis (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7663"
        }, {
          "name" : "flavonoid di-C-glucosylation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7897"
        }, {
          "name" : "hydroxymethylpyrimidine salvage",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6910"
        }, {
          "name" : "TALDO1 deficiency: failed conversion of Fru(6)P, E4P to SH7P, GA3P",
          "databaseName" : "Reactome",
          "id" : "R-HSA-6791462"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "polymethylated myricetin biosynthesis (tomato)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7160"
        }, {
          "name" : "starch degradation V",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6737"
        }, {
          "name" : "2,4,6-trinitrophenol and 2,4-dinitrophenol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7627"
        }, {
          "name" : "Carbon fixation pathways in prokaryotes",
          "databaseName" : "KEGG",
          "id" : "00720"
        }, {
          "name" : "D-glucosaminate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7310"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-73817"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Nicotinamide salvaging",
          "databaseName" : "Reactome",
          "id" : "R-SSC-197264"
        }, {
          "name" : "1,3-propanediol biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7385"
        }, {
          "name" : "lipoate biosynthesis and incorporation IV (yeast)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7382"
        }, {
          "name" : "methiin metabolism",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7614"
        }, {
          "name" : "glycogen degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5941"
        }, {
          "name" : "mandelate degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1501"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-RNO-389661"
        }, {
          "name" : "Ascorbate and aldarate metabolism",
          "databaseName" : "KEGG",
          "id" : "00053"
        }, {
          "name" : "Hereditary fructose intolerance",
          "databaseName" : "Reactome",
          "id" : "R-HSA-5657560"
        }, {
          "name" : "L-methionine biosynthesis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7977"
        }, {
          "name" : "Synthesis and degradation of ketone bodies",
          "databaseName" : "KEGG",
          "id" : "00072"
        }, {
          "name" : "Pyrimidine catabolism",
          "databaseName" : "Reactome",
          "id" : "R-CEL-73621"
        }, {
          "name" : "butirosin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7019"
        }, {
          "name" : "prodigiosin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7547"
        }, {
          "name" : "naringenin C-glucosylation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6602"
        }, {
          "name" : "10-methylstearate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8172"
        }, {
          "name" : "CMP-diacetamido-8-epilegionaminic acid biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7719"
        }, {
          "name" : "crotonosine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8100"
        }, {
          "name" : "thiamine diphosphate biosynthesis III (Staphylococcus)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6907"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-500753"
        }, {
          "name" : "L-pyrrolysine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6994"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-DDI-71336"
        }, {
          "name" : "Pentose phosphate pathway",
          "databaseName" : "Reactome",
          "id" : "R-RNO-71336"
        }, {
          "name" : "virginiae butanolide type gamma-butyrolactones biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7881"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "8-amino-7-oxononanoate biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6519"
        }, {
          "name" : "cellulose and hemicellulose degradation (cellulolosome)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6784"
        }, {
          "name" : "palmitate biosynthesis (type I fatty acid synthase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5994"
        }, {
          "name" : "UMP biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7790"
        }, {
          "name" : "dhurrin degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5976"
        }, {
          "name" : "D-xylose degradation IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7294"
        }, {
          "name" : "Peroxisomal lipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-MMU-390918"
        }, {
          "name" : "Indole alkaloid biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00901"
        }, {
          "name" : "linustatin bioactivation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7091"
        }, {
          "name" : "heme degradation V",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7846"
        }, {
          "name" : "carbaryl degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8111"
        }, {
          "name" : "Synthesis of Ketone Bodies",
          "databaseName" : "Reactome",
          "id" : "R-HSA-77111"
        }, {
          "name" : "gossypol biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5773"
        }, {
          "name" : "stearate biosynthesis II (bacteria and plants)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5989"
        }, {
          "name" : "Molybdenum cofactor biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-947581"
        }, {
          "name" : "3-hydroxy-4-methyl-anthranilate biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7765"
        }, {
          "name" : "Amaryllidacea alkaloids biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7826"
        }, {
          "name" : "taxiphyllin bioactivation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7089"
        }, {
          "name" : "5,6-dimethylbenzimidazole biosynthesis II (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7729"
        }, {
          "name" : "gliotoxin inactivation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7534"
        }, {
          "name" : "Peroxisomal protein import",
          "databaseName" : "Reactome",
          "id" : "R-RNO-9033241"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70350"
        }, {
          "name" : "Starch and sucrose metabolism",
          "databaseName" : "KEGG",
          "id" : "00500"
        }, {
          "name" : "fatty acid biosynthesis initiation (plant mitochondria)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6799"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Glycosphingolipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-SPO-1660662"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "O-Antigen nucleotide sugar biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00541"
        }, {
          "name" : "Platelet degranulation",
          "databaseName" : "Reactome",
          "id" : "R-RNO-114608"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "autoinducer AI-2 biosynthesis II (Vibrio)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6154"
        }, {
          "name" : "cuticular wax biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-282"
        }, {
          "name" : "perchlorate reduction",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6530"
        }, {
          "name" : "cis-vaccenate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5973"
        }, {
          "name" : "Pyrimidine catabolism",
          "databaseName" : "Reactome",
          "id" : "R-HSA-73621"
        }, {
          "name" : "6,7,4'-trihydroxyisoflavone biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5061"
        }, {
          "name" : "Glycosaminoglycan degradation",
          "databaseName" : "KEGG",
          "id" : "00531"
        }, {
          "name" : "L-lysine biosynthesis V",
          "databaseName" : "MetaCyc",
          "id" : "PWY-3081"
        }, {
          "name" : "sanguinarine and macarpine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5287"
        }, {
          "name" : "Hyaluronan uptake and degradation",
          "databaseName" : "Reactome",
          "id" : "R-MMU-2160916"
        }, {
          "name" : "Fructose catabolism",
          "databaseName" : "Reactome",
          "id" : "R-DME-70350"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "L-homophenylalanine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7275"
        }, {
          "name" : "cholesterol degradation to androstenedione III (anaerobic)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8151"
        }, {
          "name" : "2-hydroxypenta-2,4-dienoate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5162"
        }, {
          "name" : "cichoriin interconversion",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7057"
        }, {
          "name" : "3-dehydroquinate biosynthesis I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6164"
        }, {
          "name" : "Synthesis of Ketone Bodies",
          "databaseName" : "Reactome",
          "id" : "R-MMU-77111"
        }, {
          "name" : "Biotin metabolism",
          "databaseName" : "KEGG",
          "id" : "00780"
        }, {
          "name" : "Hyaluronan uptake and degradation",
          "databaseName" : "Reactome",
          "id" : "R-BTA-2160916"
        }, {
          "name" : "Glyoxylate metabolism and glycine degradation",
          "databaseName" : "Reactome",
          "id" : "R-DDI-389661"
        }, {
          "name" : "sulfoquinovose degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7446"
        }, {
          "name" : "CMP-3-deoxy-D-manno-octulosonate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1269"
        }, {
          "name" : "6-methoxypodophyllotoxin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5479"
        }, {
          "name" : "saframycin A biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7671"
        }, {
          "name" : "Peroxisomal lipid metabolism",
          "databaseName" : "Reactome",
          "id" : "R-HSA-390918"
        }, {
          "name" : "mevalonate pathway IV (archaea)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8125"
        }, {
          "name" : "Insulin effects increased synthesis of Xylulose-5-Phosphate",
          "databaseName" : "Reactome",
          "id" : "R-RNO-163754"
        }, {
          "name" : "mycolate biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWYG-321"
        }, {
          "name" : "Hyaluronan uptake and degradation",
          "databaseName" : "Reactome",
          "id" : "R-HSA-2160916"
        }, {
          "name" : "L-lyxonate degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7516"
        }, {
          "name" : "Polyketide sugar unit biosynthesis",
          "databaseName" : "KEGG",
          "id" : "00523"
        }, {
          "name" : "carotenoid cleavage",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6806"
        }, {
          "name" : "menaquinol-4 biosynthesis II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7998"
        }, {
          "name" : "acyl-[acyl-carrier protein] thioesterase pathway",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5142"
        }, {
          "name" : "nicotine biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5316"
        }, {
          "name" : "C5-Branched dibasic acid metabolism",
          "databaseName" : "KEGG",
          "id" : "00660"
        }, {
          "name" : "apicidin F biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7668"
        }, {
          "name" : "neoxaline biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7608"
        }, {
          "name" : "emodin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8166"
        }, {
          "name" : "L-glutamate biosynthesis V",
          "databaseName" : "MetaCyc",
          "id" : "PWY-4341"
        }, {
          "name" : "rifamycin B biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5984"
        }, {
          "name" : "starch degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-842"
        }, {
          "name" : "L-lysine biosynthesis III",
          "databaseName" : "MetaCyc",
          "id" : "PWY-2942"
        }, {
          "name" : "CMP-N-acetylneuraminate biosynthesis II (bacteria)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6139"
        }, {
          "name" : "Pyrimidine biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-500753"
        }, {
          "name" : "nevadensin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7298"
        }, {
          "name" : "kauralexin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6887"
        }, {
          "name" : "D-xylose degradation to ethylene glycol (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7178"
        }, {
          "name" : "anthocyanidin modification (Arabidopsis)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7450"
        }, {
          "name" : "Purine ribonucleoside monophosphate biosynthesis",
          "databaseName" : "Reactome",
          "id" : "R-BTA-73817"
        }, {
          "name" : "meleagrin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7607"
        } ]
      }
    },
    "locations" : [ {
      "start" : 2,
      "end" : 249,
      "hmmStart" : 3,
      "hmmEnd" : 247,
      "hmmLength" : 248,
      "hmmBounds" : "COMPLETE",
      "evalue" : 2.3E-108,
      "score" : 363.1,
      "envelopeStart" : 2,
      "envelopeEnd" : 249,
      "postProcessed" : true,
      "location-fragments" : [ {
        "start" : 2,
        "end" : 249,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 2.1E-108,
    "score" : 363.3,
    "model-ac" : "1r2rB00"
  }, {
    "signature" : {
      "accession" : "SSF51351",
      "name" : "Triosephosphate isomerase (TIM)",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "SUPERFAMILY",
        "version" : "1.75"
      },
      "entry" : {
        "accession" : "IPR035990",
        "name" : "TIM_sf",
        "description" : "Triosephosphate isomerase superfamily",
        "type" : "HOMOLOGOUS_SUPERFAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        } ]
      }
    },
    "locations" : [ {
      "start" : 3,
      "end" : 247,
      "hmmLength" : 257,
      "location-fragments" : [ {
        "start" : 3,
        "end" : 247,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 9.69E-96,
    "model-ac" : "0048198"
  }, {
    "signature" : {
      "accession" : "PF00121",
      "name" : "TIM",
      "description" : "Triosephosphate isomerase",
      "signatureLibraryRelease" : {
        "library" : "PFAM",
        "version" : "33.1"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 8,
      "end" : 245,
      "hmmStart" : 2,
      "hmmEnd" : 243,
      "hmmLength" : 243,
      "hmmBounds" : "C_TERMINAL_COMPLETE",
      "evalue" : 1.8E-86,
      "score" : 289.4,
      "envelopeStart" : 7,
      "envelopeEnd" : 245,
      "postProcessed" : true,
      "location-fragments" : [ {
        "start" : 8,
        "end" : 245,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 1.6E-86,
    "score" : 289.6,
    "model-ac" : "PF00121"
  }, {
    "signature" : {
      "accession" : "TIGR00419",
      "name" : "TIGR00419",
      "description" : "tim: triose-phosphate isomerase",
      "signatureLibraryRelease" : {
        "library" : "TIGRFAM",
        "version" : "15.0"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 8,
      "end" : 240,
      "hmmStart" : 1,
      "hmmEnd" : 228,
      "hmmLength" : 228,
      "hmmBounds" : "COMPLETE",
      "evalue" : 2.2E-77,
      "score" : 257.9,
      "envelopeStart" : 8,
      "envelopeEnd" : 240,
      "postProcessed" : false,
      "location-fragments" : [ {
        "start" : 8,
        "end" : 240,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 2.0E-77,
    "score" : 258.0,
    "model-ac" : "TIGR00419"
  }, {
    "signature" : {
      "accession" : "cd00311",
      "name" : "TIM",
      "description" : "TIM",
      "signatureLibraryRelease" : {
        "library" : "CDD",
        "version" : "3.18"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70263"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "D-apionate degradation II (RLP decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8090"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70263"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Glycolysis / Gluconeogenesis",
          "databaseName" : "KEGG",
          "id" : "00010"
        }, {
          "name" : "Entner-Doudoroff pathway I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8004"
        }, {
          "name" : "D-apionate degradation I (xylose isomerase family decarboxylase)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-8091"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352875"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-352882"
        }, {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SPO-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        }, {
          "name" : "Carbon fixation in photosynthetic organisms",
          "databaseName" : "KEGG",
          "id" : "00710"
        }, {
          "name" : "Propanoate metabolism",
          "databaseName" : "KEGG",
          "id" : "00640"
        }, {
          "name" : "Metabolic pathways",
          "databaseName" : "KEGG",
          "id" : "01100"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-HSA-70263"
        }, {
          "name" : "erythritol degradation I",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7789"
        }, {
          "name" : "superpathway of glucose and xylose degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6901"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DRE-70263"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-CEL-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-RNO-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-PFA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70263"
        }, {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Fructose and mannose metabolism",
          "databaseName" : "KEGG",
          "id" : "00051"
        }, {
          "name" : "L-threitol degradation",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7787"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-XTR-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70171"
        }, {
          "name" : "glycerol degradation to butanol",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7003"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 7,
      "end" : 246,
      "sites" : [ {
        "description" : "catalytic triad",
        "numLocations" : 3,
        "siteLocations" : [ {
          "start" : 166,
          "end" : 166,
          "residue" : "E"
        }, {
          "start" : 96,
          "end" : 96,
          "residue" : "H"
        }, {
          "start" : 14,
          "end" : 14,
          "residue" : "K"
        } ]
      }, {
        "description" : "dimer interface",
        "numLocations" : 13,
        "siteLocations" : [ {
          "start" : 83,
          "end" : 83,
          "residue" : "M"
        }, {
          "start" : 50,
          "end" : 50,
          "residue" : "D"
        }, {
          "start" : 47,
          "end" : 47,
          "residue" : "A"
        }, {
          "start" : 53,
          "end" : 53,
          "residue" : "R"
        }, {
          "start" : 86,
          "end" : 86,
          "residue" : "D"
        }, {
          "start" : 46,
          "end" : 46,
          "residue" : "T"
        }, {
          "start" : 99,
          "end" : 99,
          "residue" : "R"
        }, {
          "start" : 65,
          "end" : 65,
          "residue" : "Q"
        }, {
          "start" : 48,
          "end" : 48,
          "residue" : "Y"
        }, {
          "start" : 87,
          "end" : 87,
          "residue" : "C"
        }, {
          "start" : 98,
          "end" : 98,
          "residue" : "E"
        }, {
          "start" : 15,
          "end" : 15,
          "residue" : "M"
        }, {
          "start" : 12,
          "end" : 12,
          "residue" : "N"
        } ]
      }, {
        "description" : "substrate binding site",
        "numLocations" : 9,
        "siteLocations" : [ {
          "start" : 172,
          "end" : 172,
          "residue" : "G"
        }, {
          "start" : 231,
          "end" : 231,
          "residue" : "L"
        }, {
          "start" : 233,
          "end" : 233,
          "residue" : "G"
        }, {
          "start" : 166,
          "end" : 166,
          "residue" : "E"
        }, {
          "start" : 234,
          "end" : 234,
          "residue" : "G"
        }, {
          "start" : 212,
          "end" : 212,
          "residue" : "S"
        }, {
          "start" : 96,
          "end" : 96,
          "residue" : "H"
        }, {
          "start" : 14,
          "end" : 14,
          "residue" : "K"
        }, {
          "start" : 12,
          "end" : 12,
          "residue" : "N"
        } ]
      } ],
      "evalue" : 1.8089E-130,
      "score" : 366.476,
      "location-fragments" : [ {
        "start" : 7,
        "end" : 246,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "cd00311"
  }, {
    "signature" : {
      "accession" : "PTHR21139:SF24",
      "name" : "TRIOSEPHOSPHATE ISOMERASE",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "PANTHER",
        "version" : "15.0"
      },
      "entry" : null
    },
    "locations" : [ {
      "start" : 1,
      "end" : 249,
      "hmmStart" : 1,
      "hmmEnd" : 249,
      "hmmLength" : 249,
      "hmmBounds" : "COMPLETE",
      "envelopeStart" : 1,
      "envelopeEnd" : 249,
      "location-fragments" : [ {
        "start" : 1,
        "end" : 249,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 0.0,
    "familyName" : "Not available",
    "score" : 585.3,
    "model-ac" : "PTHR21139:SF24"
  } ],
  "xref" : [ {
    "name" : "CAA49379.1 triosephosphate isomerase [Homo sapiens]",
    "id" : "CAA49379.1"
  } ]
} ]
}'''


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
    testpro = '''>QAB35645.1 BRI1 protein [Solanum tuberosum]
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
    ips.application_select(['TIGRFAM', 'CDD', 'PfamA'])
    # ips.application_select(['Phobius', 'ProSitePatterns'])
    ips.output_select('json')
    ips.parameter_select({'goterms': True, 'pathways': True})

    if not ips.run():
        exit(1)

    time.sleep(ips.poll_time)
    while ips.status() != 'FINISHED':
        ips.poll_count += 1
        if ips.poll_count > ips.poll_max:
            break

        time.sleep(ips.poll_time)

    ips.result()

    # parse and print the result
    # ips.content = json_test()
    parsed_result = ips.parse_json()

    for eachmotif in parsed_result['motifs']:
        print('{ipr_accession}\t{src_accession}\t{description}'.format(**eachmotif))
    for goterm in parsed_result['go']:
        go = parsed_result['go'][goterm]
        print('{}\t{}\t{}'.format(goterm, go['name'], go['source']))
    for path in parsed_result['pathway']:
        pathway = parsed_result['pathway'][path]
        print('{}\t{}\t{}'.format(path, pathway['name'], pathway['source']))

    print('done')

    exit(0)

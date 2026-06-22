import sys
import time
import json
import requests
from api.jobmanager_api import JobManagerAPI


class InterproscanAPI(JobManagerAPI):
    """=============================================================================================
    Container for the methods that must be provided to JobManagerAPI
    ============================================================================================="""

    def submit(self, show_query=False):
        """-----------------------------------------------------------------------------------------
        Construct a REST command and dispatch the job to the server
        Any previously existing jobID is overwritten

        :param show_query: boolean, print query if true
        :return: logical, True = success, False = failure
        -----------------------------------------------------------------------------------------"""
        is_success = False

        # general fields for all queries
        param = {u'program': 'iprscan6', u'email': self.email, u'title': self.title, u'sequence': self.sequence,
                 u'output' : self.output, u'stype': 'p'}

        if self.applications:
            # add selected applications
            param['appl'] = ','.join(self.applications)

        if self.parameters:
            # add the goterms and pathways parameters
            for para in self.parameters:
                param[para] = self.parameters[para]

        command = self.url + 'run'
        self.response = requests.post(command, files=param, headers={'User-Agent': 'ips-client'})

        if show_query:
            # print out query if requested
            print(self.response.request.headers, '\n')
            print(self.response.request.body, '\n')

        if self.response_is_error('submitting job'):
            self.jobstatus = 'failed'
        else:
            # success
            self.jobid = self.response.text
            self.jobstatus = 'submitted'
            self.message = {'type'    : 'submitted',
                            'text'    : f'job_name={self.title};job_id={self.jobid}',
                            'loglevel': 1}

            is_success = True

        return is_success

    def status(self, log=True):
        """-----------------------------------------------------------------------------------------
        Poll job status at the server. The job is polled only once so if you want to poll
        multiple times call this method in a loop

        :return: string, status of job at server
        -----------------------------------------------------------------------------------------"""
        command = self.url + 'status/' + self.jobid
        self.response = requests.get(command)
        response_text = self.response.text.rstrip()

        if 'RUNNING' in self.response.text:
            self.jobstatus = 'running'
            self.message = {'type'    : 'polling',
                            'text'    : f'job_id={self.jobid};response={response_text}',
                            'loglevel': 2}

        elif 'FINISHED' in self.response.text:
            if self.jobstatus != 'finished':
                # only print finished message once
                self.jobstatus = 'finished'
                self.message = {'type'    : 'finished',
                                'text'    : f'job_id={self.jobid}',
                                'loglevel': 1}

        return self.jobstatus

    def result(self):
        """-----------------------------------------------------------------------------------------
        Retrieve the result

        :return: string, 'retrieved' if successful, '' if unsuccessful (False)
        -----------------------------------------------------------------------------------------"""
        # get the final result
        command = self.url + 'result/' + self.jobid + '/' + self.output
        self.response = requests.get(command)
        if not self.response_is_error('retrieving result'):
            # success
            self.content = self.response.text
            self.message = {'type'    : 'retrieved',
                            'text'    : f'job_id={self.jobid};output_len={len(self.output)}',
                            'loglevel': 1}
            return 'retrieved'

        return ''


class InterproscanQuery():
    """=============================================================================================
    Interpro class for running interproscan.

    This class holds the information needed to run interproscan via the iprscan6 API. Management of
    job quest is in the JobManagerAPI. The Interpro object corresponds to one query.

    25 December 2018    Michael Gribskov
    ============================================================================================="""
    # available options/parameters taken from
    # https://www.ebi.ac.uk/jdispatcher/docs/webservices/

    # available outputs
    # from https://www.ebi.ac.uk/Tools/services/rest/iprscan5/resulttypes
    # log - The output from the tool itself
    # out - The results of the job (XML format)
    # tsv - The results of the job in text format, tab separated values
    # xml - The results of the job in XML
    # gff - The results of the job in GFF3 format (gff3 is a synonym)
    # json - The results of the job in JSON format
    # jsonl - JSON format intended for streaming
    # sequence - Input sequence as seen by the tool
    # submission - The submission details which were submitted as a job
    # zip - full result zipped
    available = {'applications': ['AntiFam', 'CATH-Gene3D', 'CATH-FunFam', 'CDD', 'COILS', 'HAMAP',
                                  'MobiDB-lite', 'NCBIFAM', 'PANTHER', 'Pfam', 'Phobius', 'PIRSF',
                                  'PRINTS', 'PROSITE-patterns', 'PROSITE-profiles', 'SFLD', 'SMART',
                                  'SUPERFAMILY', 'SignalP-Euk', 'SignalP-Prok'],
                 'commands'    : ['run', 'status', 'result'],
                 'outputs'     : ['out', 'log', 'tsv', 'xml', 'gff', 'gff3', 'error', 'json', 'jsonl',
                                  'submission', 'zip']
                 }

    def __init__(self, log_level=1):
        """-----------------------------------------------------------------------------------------
        interpro query/response constructor

        loglevel   0 no log, 1 job submission/completion, 2 all
        -----------------------------------------------------------------------------------------"""

        self.email = ''  # user email (optional)
        self.title = ''  # title for job (optional)
        self.sequence = ''
        self.applications = []
        self.output = 'json'
        self.parameters = {}
        self.log_level = log_level

        self.url = u'https://www.ebi.ac.uk/Tools/services/rest/iprscan6/'
        self.jobid = ''
        self.jobname = ''
        self.jobstatus = ''

        self.response = None
        self.content = ''

    def application_select(self, selected, keep=False):
        """-----------------------------------------------------------------------------------------
        Add a list of applications to be run.  Each application in the list is compared to the
        available applications and if not present a warning is issued.  The default is to run all
        applications so an empty selected list signifies the default.

        :param selected: list of strings, selected applications to run
        :param keep: Boolean, retain current applications, just add new ones
        :return: int, number of selected applications
        -----------------------------------------------------------------------------------------"""
        if not keep:
            self.applications = []

        for app in selected:
            # if app == 'Pfam':
            #     app = 'PfamA'
            if app in InterproscanQuery.available['applications']:
                self.applications.append(app)
            else:
                self.message = {'type'    : 'not_available',
                                'text'    : f'application={app}',
                                'loglevel': 2}

        return len(self.applications)

    def output_select(self, selected):
        """-----------------------------------------------------------------------------------------
        select the output format.  Only one is allowed

        :param selected: string, one of the formats in self.output_avail
        :return: True if format is available
        -----------------------------------------------------------------------------------------"""
        # self.output = ''
        if selected in InterproscanQuery.available['outputs']:
            self.output = selected
        else:
            self.message = {'type'    : 'not_available',
                            'text'    : f'output={selected}',
                            'loglevel': 2}

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
                      'category': ontology category = BIOLOGICAL_PROCESS |
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
                           'description'  : entry['description'] or ''})

            # name = entry['name']
            # type = entry['type']

            if 'goXRefs' in entry:
                gostr = ''
                for go in entry['goXRefs']:
                    gostr += '{} ({}:{})'.format(go['id'], go['category'], go['name'])
                    if go['id'] in go_all:
                        go_all[go['id']]['source'].append(source_accession)
                    else:
                        go_all[go['id']] = {'name'  : go['name'], 'category': go['category'],
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
                        path_all[id] = {'name'  : path['name'],
                                        'source': [source_accession]}

        return {'jobname': self.jobname, 'motifs': motifs, 'go': go_all, 'pathway': path_all}

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
            self.jobstatus = 'error'
            self.message = {'type'    : task,
                            'text'    : f'job_id={self.jobid};status={self.response.status_code}',
                            'loglevel': 1}

        return is_error


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    ips = InterproscanQuery()
    query = {}
    ips.email = 'gribskov@purdue.edu'
    # ips.title = 'BRI1'
    # ips.sequence = testseq[1]
    ips.application_select(['AntiFam', 'SignalP-Euk', 'Pfam'])
    ips.output_select('gff')
    ips.parameter_select({'goterms': True, 'pathways': False})

    print('submitting')
    if not ips.submit():
        exit(1)

    poll_time = 20
    poll_count = 0
    poll_max = 50
    while ips.status() != 'finished':
        time.sleep(poll_time)
        poll_count += 1
        print(f'\t ... polling({poll_count}) = {ips.jobstatus}')
        if poll_count > poll_max:
            break

    print('collecting result')
    ips.result()

    # parse and print the result, comment out the above and ncomment the next line to test parsing
    # the result without running a query
    # ips.content = json_test()
    if ips.output.startswith('json'):
        parsed_result = ips.parse_json()

        for eachmotif in parsed_result['motifs']:
            print('{ipr_accession}\t{src_accession}\t{description}'.format(**eachmotif))
        for goterm in parsed_result['go']:
            go = parsed_result['go'][goterm]
            print('{}\t{}\t{}'.format(goterm, go['name'], go['source']))
        for path in parsed_result['pathway']:
            pathway = parsed_result['pathway'][path]
            print('{}\t{}\t{}'.format(path, pathway['name'], pathway['source']))
    else:
        print(f'{ips.content}')

    print('done')

    exit(0)

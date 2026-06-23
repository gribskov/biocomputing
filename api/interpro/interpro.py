import sys
import time
import json
import requests
from api.jobmanager_api import JobManagerAPI


class InterproscanAPI(JobManagerAPI):
    """=============================================================================================
    Container for the methods that must be provided to JobManagerAPI
    ============================================================================================="""

    def submit(self, query, show_query=False):
        """-----------------------------------------------------------------------------------------
        Construct a REST command and dispatch the job to the server
        Any previously existing jobID is overwritten

        :param query: InterproscanQuery     query object
        :param show_query: boolean          print query if true
        :return: logical                    True = success, False = failure
        -----------------------------------------------------------------------------------------"""
        is_success = False

        query_fields = ['email', 'title', 'goterms', 'pathways', 'stype', 'appl', 'sequence']
        param = {}
        for q in query_fields:
            param[q] = query.parameters[q]

        command = query.parameters['url'] + 'run'
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
            self.message = {'type': 'submitted',
                            'text': f'job_name={self.title};job_id={self.jobid}',
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
            self.message = {'type': 'polling',
                            'text': f'job_id={self.jobid};response={response_text}',
                            'loglevel': 2}

        elif 'FINISHED' in self.response.text:
            if self.jobstatus != 'finished':
                # only print finished message once
                self.jobstatus = 'finished'
                self.message = {'type': 'finished',
                                'text': f'job_id={self.jobid}',
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
            self.message = {'type': 'retrieved',
                            'text': f'job_id={self.jobid};output_len={len(self.output)}',
                            'loglevel': 1}
            return 'retrieved'

        return ''

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
            self.message = {'type': task,
                            'text': f'job_id={self.jobid};status={self.response.status_code}',
                            'loglevel': 1}

        return is_error


class InterproscanQuery():
    """=============================================================================================
    Interpro class for running interproscan.

    This class holds the information needed to run interproscan via the iprscan6 API. Management of
    job quest is in the JobManagerAPI. The Interpro object corresponds to one query.

    available options/parameters taken from https://www.ebi.ac.uk/jdispatcher/docs/webservices/

    available outputs
    from https://www.ebi.ac.uk/Tools/services/rest/iprscan5/resulttypes
    log - The output from the tool itself
    out - The results of the job (XML format)
    tsv - The results of the job in text format, tab separated values
    xml - The results of the job in XML
    gff - The results of the job in GFF3 format (gff3 is a synonym)
    json - The results of the job in JSON format
    jsonl - JSON format intended for streaming
    sequence - Input sequence as seen by the tool
    submission - The submission details which were submitted as a job
    zip - full result zipped

    25 December 2018    Michael Gribskov
    ============================================================================================="""
    # the keywords agree with the Interpro definition of query fields
    available = {'appl': ['AntiFam', 'CATH-Gene3D', 'CATH-FunFam', 'CDD', 'COILS', 'HAMAP',
                          'MobiDB-lite', 'NCBIFAM', 'PANTHER', 'Pfam', 'Phobius', 'PIRSF',
                          'PRINTS', 'PROSITE-patterns', 'PROSITE-profiles', 'SFLD', 'SMART',
                          'SUPERFAMILY', 'SignalP-Euk', 'SignalP-Prok'],
                 'command': ['run', 'status', 'result'],
                 'resultType': ['out', 'log', 'tsv', 'xml', 'gff', 'gff3', 'error', 'json', 'jsonl',
                                'submission', 'zip']
                 }

    def __init__(self, query_fields):
        """ - --------------------------------------------------------------------------------------
        variables needed for interproscan query
        email
        title
        sequence
        stype
        output

        loglevel
            0 no log,
            1 job submission / completion
            2 all
        -----------------------------------------------------------------------------------------"""
        self.parameters = query_fields

        self.jobid = ''
        self.jobstatus = ''
        self.response = ''
        self.content = ''

    def validate(self, keys):
        """-----------------------------------------------------------------------------------------
        Validate query information vs the list of available options. Only implemented for
        'applications' and 'output'

        :param keys: list       list of option keys to validate (keys in available)
        :return: str            description of actions
        -----------------------------------------------------------------------------------------"""
        errstr = ''
        resulttype_default = 'gff3'
        for key in keys:
            if key == 'appl':
                valid_appl = []
                for k in self.parameters['appl']:
                    if k in self.available['appl']:
                        valid_appl.append(k)
                    else:
                        errstr += f'application {k} invalid\n'
                self.parameters['appl'] = valid_appl

            elif key == 'resultType':
                if self.parameters['resultType'] not in InterproscanQuery.available['resultType']:
                    self.parameters['resultType'] = resulttype_default
                    errstr += f'resultType {self.parameters['resultType']} not supported, '
                    errstr += f'resultType set to {resulttype_default}\n'

        return errstr

    def parse_json(self):
        """-----------------------------------------------------------------------------------------
        Parse the content returned from the server in JSON format. Three outputs are produced
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
                           'description': entry['description'] or ''})

            # name = entry['name']
            # type = entry['type']

            if 'goXRefs' in entry:
                gostr = ''
                for go in entry['goXRefs']:
                    gostr += '{} ({}:{})'.format(go['id'], go['category'], go['name'])
                    if go['id'] in go_all:
                        go_all[go['id']]['source'].append(source_accession)
                    else:
                        go_all[go['id']] = {'name': go['name'], 'category': go['category'],
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
                        path_all[id] = {'name': path['name'],
                                        'source': [source_accession]}

        return {'jobname': self.jobname, 'motifs': motifs, 'go': go_all, 'pathway': path_all}




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

class Interpro:
    """=============================================================================================
    Interpro class for running interproscan

    25 December 2018    Michael Gribskov
    ============================================================================================="""
    import sys
    import requests

    def __init__(self, loglevel=0, poll_time=60, poll_count=20):
        """-----------------------------------------------------------------------------------------
        interpro query/response constructor

        track   0 no log, 1 job submission/completion, 2 all
        -----------------------------------------------------------------------------------------"""
        self.log = loglevel
        self.commands_avail = {'run', 'status', 'result'}
        self.output_avail = {'xml', 'json', 'tsv', 'out', 'gff', 'svg'}
        self.output = 'tsv'
        self.poll_time = poll_time  # seconds between polling
        self.poll_count = poll_count  # maximum number of times to poll

        self.email = ''  # user email (optional)
        self.title = ''  # title for job (optional)
        self.sequence = ''

        self.url = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/'
        self.jobid = ''
        self.state = 'UNKNOWN'
        self.result = ''

    def run(self):
        """-----------------------------------------------------------------------------------------
        Construct a REST command and dispatch the job to the server
        Any previously existing jobID is overwritten
        :return:
        -----------------------------------------------------------------------------------------"""
        # send the initial query
        param = {'email': self.email, 'title': self.title, 'sequence': self.sequence}
        command = self.url + 'run'
        response = requests.post(command, param)
        self.jobid = response.text
        # TODO add error trapping

        if self.log:
            # TODO add time and sequence name
            sys.stderr.write('job {} submitted\n'.format(self.jobid))

        return jobid

    def status(self):
        """-----------------------------------------------------------------------------------------
        Poll job status at the server

        :return: Logical, True if a result was returned
        -----------------------------------------------------------------------------------------"""
        complete = False
        tries = 0
        while not complete:
            # don't poll too often
            time.sleep(self.poll_time)
            tries += 1

            command = self.url + self.jobid
            response = requests.get(command)
            if self.log > 1:
                # TODO add time and jobid
                sys.stderr.write('    polling... response->{}\n'.format(response.text))

            if 'FINISHED' in response.text:
                complete = True
                break
            elif tries >= self.poll_count:
                break

        self.state = response.text
        if not complete:
            # polling reached limit
            # TODO add time
            if self.log > 0:
                sys.stderr.write(
                    ('unable to find result () in {} tries\n'.format(self.jobid, self.poll_count))

                else:
                # TODO add time
                if self.log > 0:
                    sys.stderr.write('interproscan {} finished\n'.format(self.jobid))

        return complete

    def result(self):
        """-----------------------------------------------------------------------------------------
        Retrieve the result

        :return:
        -----------------------------------------------------------------------------------------"""
        # get the final result
        command = '{}{}/{}/{}'.format(self.url, 'result', self.jobid, self.output)
        response = requests.get(command)
        self.result = response.text
        # TODO error trapping

        return


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    ips = Interpro()
    ips.email = 'gribskov@purdue.edu'
    ipos.title = 'globin'
    ips.sequence = '''>sp|P69905.2|HBA_HUMAN RecName: Full=Hemoglobin subunit alpha; AltName: Full=Alpha-globin; AltName: Full=Hemoglobin alpha chain
    MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNA
    VAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSK
    YR'''

    ips.run()
    if ips.status:
        ips.result()

exit(0)

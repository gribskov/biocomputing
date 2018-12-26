import sys
import time
import requests


class Interpro:
    """=============================================================================================
    Interpro class for running interproscan

    25 December 2018    Michael Gribskov
    ============================================================================================="""

    def __init__(self, loglevel=0, poll_time=60, poll_count=20):
        """-----------------------------------------------------------------------------------------
        interpro query/response constructor

        track   0 no log, 1 job submission/completion, 2 all
        -----------------------------------------------------------------------------------------"""
        self.log = loglevel
        self.commands_avail = ['run', 'status', 'result']
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
        self.content = ''

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
        if response.status_code == 200:
            # success
            self.jobid = response.text
        else:
        # failure

        #  TODO add error trapping

        if self.log:
            # TODO add sequence name?
            sys.stderr.write(
                '{}\tinterproscan job {} submitted\n'.format(Interpro.logtime(), self.jobid))

        return self.jobid

    def status(self):
        """-----------------------------------------------------------------------------------------
        Poll job status at the server

        :return: Logical, True if a result was returned
        -----------------------------------------------------------------------------------------"""
        response = None
        complete = False
        tries = 0
        while not complete:
            # don't poll too often
            time.sleep(self.poll_time)
            tries += 1

            command = self.url + 'status/' + self.jobid
            response = requests.get(command)
            if self.log > 1:
                sys.stderr.write('{}\tinterproscan job{} polling - {}\n'.format(
                    Interpro.logtime(), self.jobid, response.text))

            if 'FINISHED' in response.text:
                complete = True
                break
            elif tries >= self.poll_count:
                break

        self.state = response.text
        if not complete:
            # polling reached limit
            if self.log > 0:
                sys.stderr.write('{}\tinterproscan job {} unsuccessful after {} tries\n'.format(
                    Interpro.logtime(), self.jobid, self.poll_count))

        else:
            if self.log > 0:
                sys.stderr.write(
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
        if response.status_code == 200:
            # success
            self.content = response.text
            if self.log > 1:
                sys.stderr.write('{}\tinterproscan job {} result retrieved from {} as {}'.format(
                    Interpro.logtime(), self.jobid, self.url, self.output))
            return True

        else:
            # failure
            sys.stderr.write('{}\t{}interproscan job {} error retrieving result\tstatus={}'.format(
                Interpro.logtime(), self.jobid, response.status_code))

        return

    @classmethod
    def logtime(cls):
        """-----------------------------------------------------------------------------------------
        Return current time as a string. Format is 10/Oct/2000:13:55:36 which is similar to the
        common log format (without the time zone)

        :return: string
        -----------------------------------------------------------------------------------------"""
        return time.strftime('%d/%b/%G:%H:%M:%S', time.localtime(time.time()))

    def response_is_error(self, task, response):
        """-----------------------------------------------------------------------------------------
        Return true if the response code is other than 200. Write error message to stderr if
        loglevel > 1. Task is a string describing the task that failed for inclusion in the error message

        :param response: requests object response
        :return: logical True = error, False = no error
        -----------------------------------------------------------------------------------------"""
        is_error = False
        if response.status_code == 200:
        # success
        else:
            if self.log > 0:
                sys.stderr.write('{}\t{}interproscan job {} error {}\tstatus={}'.format(
                    Interpro.logtime(), self.jobid, task, response.status_code))
            is_error = True

        return is_error


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    ips = Interpro(loglevel=2, poll_time=20)
    ips.email = 'gribskov@purdue.edu'
    ips.title = 'globin'
    ips.sequence = '''>sp|P69905.2|HBA_HUMAN RecName: Full=Hemoglobin subunit alpha; AltName: Full=Alpha-globin; AltName: Full=Hemoglobin alpha chain
    MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNA
    VAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSK
    YR'''

    ips.run()
    if ips.status():
        ips.result()

exit(0)

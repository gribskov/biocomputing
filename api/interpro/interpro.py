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
        self.log_fh = sys.stderr
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
        :return: logical, True = success, False = failure
        -----------------------------------------------------------------------------------------"""
        is_success = False

        # send the initial query
        param = {'email': self.email, 'title': self.title, 'sequence': self.sequence}
        command = self.url + 'run'
        response = requests.post(command, param)
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
    ips = Interpro(loglevel=2, poll_time=20)
    ips.email = 'gribskov@purdue.edu'
    ips.title = 'globin'
    ips.sequence = '''>leghemoglobin [Medicago sativa]
MQIQIAKQKQKNKKRNMGFTEKQEALVNSSFESFKQNPGYSVLFYTIILEKAPAAKGMFSFLKDSAGVQD
SPKLQAHAGKVFGMVRDSAAQLRATGGVVLGDATLGAIHIQNGVVDPHFVVVKEALLKTIKESSGDKWSE
ELSTAWEVAYDALATAIKKAMS
'''

    if not ips.run():
        exit(1)

    if ips.status():
        ips.result()
        sys.stdout.write(ips.content)

    exit(0)

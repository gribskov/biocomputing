import sys
import time
from api.blast.blastncbi import BlastNCBI
from api.interpro.interpro import Interpro


class Jobmanager:
    """=============================================================================================
    jobmanager

    manage submission of jobs using defined APIs. APIs should implement
        submit
        status
        result
        clone
        poke

    Michael Gribskov     02 April 2021
    ============================================================================================="""

    def __init__(self, loglevel=1, log_fh=sys.stderr, poll_time=61, poll_wait=0, poll_max=50):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""

        self.joblist = []
        self.api = None
        self.template = None

        self.loglevel = loglevel
        self.log_fh = log_fh
        self.poll_time = poll_time  # seconds between polling cycles
        self.poll_wait = poll_wait  # seconds between polling queries
        self.poll_max = poll_max  # maximum number of times to poll
        self.poll_count = 0  # number of times this job has been polled

    def log_message(self, type, message, loglevel=1):
        """-----------------------------------------------------------------------------------------
        write a message to the log. Messages are only written if self.loglevel >= loglevel

        types:
            not_available (for submission options)
            submitted
            polling
            finished
            retrieved
            server_error

        :type: string, type of message
        :param message: string, text of message
        :param loglevel: int, minimum loglevel to write message to log
        :return: True if message was written
        -----------------------------------------------------------------------------------------"""
        if self.loglevel >= loglevel:
            event_time = time.strftime('%d/%b/%G:%H:%M:%S', time.localtime(time.time()))
            self.log_fh.write('{}\t{}\t{}\n'.format(event_time, type, message))
            return True

        return False

    def new_job(self, *args, **kwds):
        """-----------------------------------------------------------------------------------------
        Return a new object from the current API.  pass positional and keyword args through to
        the object constructor

        :return: object from self.api
        -----------------------------------------------------------------------------------------"""
        return self.api(*args, **kwds)

    def new_job_from_template(self):
        """-----------------------------------------------------------------------------------------
        Return a copy of the template object.  The object must implement clone().
        If template has not been loaded, return the default object

        :return: object from self.template
        -----------------------------------------------------------------------------------------"""
        if self.template:
            return self.template.clone()
        else:
            return self.api()

    def poll(self, joblist, wait_all=True):
        """---------------------------------------------------------------------------------------------
        Poll all jobs in the jobs list and record their status. If wait_all is True, poll until all
        have finished.

        :param wait_all: boolean, if True poll until all jobs finish
        :return: int, number of finished jobs
        ---------------------------------------------------------------------------------------------"""
        time.sleep(self.poll_time)

        not_done = True
        n = 0
        n_finished = 0
        while not_done:
            not_done = False
            n += 1

            for job in self.joblist:
                self.joblist[job] = job.status()
                if job.status() == 'finishe':
                    n_finished += 1
                time.sleep(self.poll_wait)

            if wait_all:
                not_done = n != n_finished

        return n_finished

    def save_all(self, reformat=None, fh=None, remove=True):
        """---------------------------------------------------------------------------------------------
        Return the output of all finished jobs.
        Reformat is a callback function used to reformat the output.  for instance, interpro.parse_json
        If fh is True, output is written to the filehandle after reformatting.
        If remove is true, jobs are deleted from the list after saving

        :param joblist: dict, ips object is key, staus is value
        :param reformat: function, callback function for formatting job result, argument is ips object
        :param fh: filehandle for writable file
        :param remove: boolean, remove finished jobs after saving
        :return: string, text of job content
        ---------------------------------------------------------------------------------------------"""
        delete_list = []
        for job in joblist:
            if joblist[job] != 'finished':
                # skip unfinished jobs
                continue

            joblist[job] = job.result()  # retrieve the completed job
            # text = job.content
            if reformat:
                text = reformat(job)

            if fh:
                if text:
                    fh.write('!{} - {}s\n'.format(job.jobname, job.jobid))
                    fh.write('{}\n'.format(text))
                else:
                    fh.write('!{} - {} no hits\n'.format(job.jobname, job.jobid))

            if remove:
                delete_list.append(job)

        for job in delete_list:
            del joblist[job]

        return text

    def submit(self, job):
        """-----------------------------------------------------------------------------------------
        Submit a job through the self.api and add the job to the joblist.  The object is the key and
        the value is the status

        :param job: object from api
        :return: int, number of jobs in the joblist
        -----------------------------------------------------------------------------------------"""
        status = job.submit()
        self.joblist[job] = status

        self.log_message(status, 'job_name={};job_id={}'.format(job.jobname, job.jobid), loglevel=1)

        return len(self.joblist)


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    joblist = Jobmanager()
    joblist.api = BlastNCBI
    print(joblist.api.poke())

    new = joblist.api()

    exit(0)

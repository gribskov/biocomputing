from abc import ABC, abstractmethod
import time
import sys


class JobManagerAPI(ABC):
    """=============================================================================================
    Abstract base class for jobmanager API

    Children inherit
        clone
        poke

    Children must implement
        submit - returns self.rid
        status - returns self.jobstatus
        result - returns the final result

    children inherit attributes
        message - transmit log message to manager
        jobname - name of job for use as key in final parsed result
        content - retrieved content of result

    Michael Gribskov     19 April 2021
    ============================================================================================="""
    # class variables shared between all instances
    poll_delay = 60
    poll_maxcount = 25
    simultaneous_jobs = 1
    joblist = []

    # ----------------------------------------------------------------------------------------------
    # concrete methods available to all subclasses
    # ----------------------------------------------------------------------------------------------
    def __init_subclass__(cls):
        """-----------------------------------------------------------------------------------------
        called when subclass is instantiated.
        message allows the subclass to send error/log output to the manager.  Allowed types
            not_available (for submission options)
            failed
            submitted
            polling
            finished
            retrieved
            server_error

            joblist: dict   key is a unique id, value is status
        :return:
        -----------------------------------------------------------------------------------------"""
        cls.jobname = ''
        cls.message = []
        cls.content = ''
        # cls.poll_delay = 10
        # cls.poll_maxcount = 25
        # cls.simultaneous_jobs = 1
        # cls.joblist = []
        pass

    def clone(self):
        """-----------------------------------------------------------------------------------------
        Return a copy of the object.  This allows an object with the metadata filled in to be
        used as a template for a series of jobs

        INHERITABLE

        :return: copy of object
        -----------------------------------------------------------------------------------------"""
        # type(self)() correctly targets the subclass (e.g., Interpro)
        copy = type(self)()
        for v in vars(self):
            setattr(copy, v, getattr(self, v))

        return copy

    def poke(self):
        """-----------------------------------------------------------------------------------------
        Return a signature string.  Useful to identify the class when it is used as a callback

        INHERITABLE

        :return: string
        -----------------------------------------------------------------------------------------"""
        return self.__class__.__name__

    # ----------------------------------------------------------------------------------------------
    # abstract methods to be supplied by the subclass
    # ----------------------------------------------------------------------------------------------

    @abstractmethod
    def result(self, *args, **kwds):
        """-----------------------------------------------------------------------------------------
        Retrieve the result

        -----------------------------------------------------------------------------------------"""

    @abstractmethod
    def status(self, *args, **kwds):
        """-----------------------------------------------------------------------------------------
        Checks to see if job is complete.  Often this means polling the server and getting the
        status of the job in self.rid.

        self.status is set as running, unknown, or finished

        :return: string     status: running | unknown | finished
        -----------------------------------------------------------------------------------------"""

    @abstractmethod
    def submit(self, *args, **kwds):
        """-----------------------------------------------------------------------------------------
        Start the job, usually by submitting to the service
        :return: string, request ID (rid)
        -----------------------------------------------------------------------------------------"""

    # -----------------------------------------------------------------------------------------------
    # high level methods built on the methods supplied by the subclass
    # ----------------------------------------------------------------------------------------------

    def poll_all(self):
        """-----------------------------------------------------------------------------------------
        Poll the jobs in the jobs_pending list repeatedly until all have finished. Finished can be
            1) reached maximum number of polling attempts
            2) returned a status other than success or waiting
            3) success

        uses class variables:
            joblist: list of interproscan objects that have been submitted
            poll_delay: int, seconds to wait between polling
            poll_maxcount: int, maximum number of times to poll
            ntries: int, number of polling trials
        -----------------------------------------------------------------------------------------"""
        joblist = self.joblist
        ntries = 0
        jobs_running = True
        while jobs_running:
            ntries += 1
            jobs_running = False
            time.sleep(self.poll_delay)

            # the jobs in joblist are InterproscanQuery objects
            for job in joblist:
                if job.jobstatus != 'finished':
                    status = self.status(job)
                    # job.jobstatus = status
                    if status == 'running':
                        jobs_running = True
                        # no need to keep checking after one running job is found
                        # return to while jobs_running
                        break

        return ntries

    def result_all(self):
        """-----------------------------------------------------------------------------------------
        Retrieve results for all queries marked as finished using the result() method provided by
        the concrete subclass.

        uses class variables:
            joblist: list of interproscan objects that have been submitted
        uses subclass result() method
        -----------------------------------------------------------------------------------------"""
        joblist = self.joblist
        # the jobs in joblist are InterproscanQuery objects
        for job in joblist:
            if job.jobstatus == 'finished':
                status = self.result(job)

        return status

    def save_all(self):
        """-----------------------------------------------------------------------------------------
        Retrieve results for all queries marked as finished. Output file name is based on the query
        jobname

        uses class variables:
            joblist: list of interproscan objects that have been submitted
        uses subclass result() method
        -----------------------------------------------------------------------------------------"""
        joblist = self.joblist
        # the jobs in joblist are InterproscanQuery objects
        nwritten = 0
        for job in joblist:
            if job.jobstatus == 'finished':
                if job.parameters['jobname'] == 'stdout':
                    outfile = sys.stdout
                else:
                    outfile = open(f'{job.parameters['jobname']}.out', 'w')

                outfile.write(job.content)
                outfile.close()
                nwritten += 1

        return nwritten


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # should succeed
    class success(JobManagerAPI):

        def __init__(self):
            self.id = 'success'

        def submit(self):
            pass

        def status(self):
            pass

        def result(self):
            pass


    # should fail, does not implement abstract methods
    class fail(JobManagerAPI):

        def __init__(self):
            self.id = 'fail'


    test = success()
    print(f'{test.id}, class={test.poke()}\n')

    test = fail()
    print(f'{test.id}, class={test.poke()}')
    test.submit()

    exit(0)

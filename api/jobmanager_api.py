from abc import ABC, abstractmethod


class JobManagerAPI(ABC):
    """=============================================================================================
    Abstract base class for jobmanager API

    Children inherit
        clone
        poke

    Children must implement
        submit - returns self.rid
        status - reeturns self.jobstatus
        result

    Michael Gribskov     19 April 2021
    ============================================================================================="""

    def clone(self):
        """-----------------------------------------------------------------------------------------
        Return a copy of the object.  This allows an object with the metadata filled in to be
        used as a template for a series of jobs

        INHERITABLE

        :return: copy of object
        -----------------------------------------------------------------------------------------"""
        copy = __class__()
        for v in vars(self):
            setattr(copy, v, getattr(self, v))

        return copy

    def poke(self):
        """-----------------------------------------------------------------------------------------
        Return a signature string.  Useful to identify the when class it is used as a callback

        INHERITABLE

        :return: string
        -----------------------------------------------------------------------------------------"""
        return __class__.__name__

    @abstractmethod
    def result(self):
        """-----------------------------------------------------------------------------------------
        Retrieve the result

        :return: None
        -----------------------------------------------------------------------------------------"""
        # return None

    @abstractmethod
    def status(self, wait=False):
        """-----------------------------------------------------------------------------------------
        Checks to see if job is complete.  Often this means polling the server and getting the
        status of the job in self.rid.

        self.status is set as running, unknown, or finished, respectively

        :param wait: boolean, if True wait for self.poll_time before sending request
        :return: string, status running | unknown | finished
        -----------------------------------------------------------------------------------------"""
        return self.jobstatus

    @abstractmethod
    def submit(self):
        """-----------------------------------------------------------------------------------------
        Start the job, usually by submitting to the service
        :return: string, request ID (rid)
        -----------------------------------------------------------------------------------------"""
        # return self.rid


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


    class fail(JobManagerAPI):

        def __init__(self):
            self.id = 'fail'


    test = success()
    print(f'{test.id}, class={test.poke()}\n')

    test = fail()
    print(f'{test.id}, class={test.poke()}')
    test.submit()

    exit(0)

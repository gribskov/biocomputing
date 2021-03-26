"""=================================================================================================
Run interproscan for multiple sequences using the EMBL interproscan server

Usage
    intropro_batch.py <Fasta file>

26 December 2018    Michael Gribskov
================================================================================================="""
import sys
import time
import argparse
from sequence.fasta import Fasta
from api.interpro.interpro import Interpro

import textwrap as _textwrap


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        self._max_help_position = 30
        return _textwrap.wrap(text, 60)


def arguments_get():
    minlen_default = 80
    loglevel_default = 1
    batch_limit_default = 20
    batch_wait_default = 60
    cl = argparse.ArgumentParser(description='Interproscan of ORF sequences',
                                 formatter_class=CustomFormatter)
    cl.add_argument('--logfile', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output file for log information')
    cl.add_argument('-m', '--minlen', type=int, default=minlen_default,
                    help='Minimum length ORF to run')
    cl.add_argument('--batch_limit', type=int, default=batch_limit_default,
                    help='Number of sequences to submit per batch')
    cl.add_argument('--batch_wait', type=int, default=batch_wait_default,
                    help='Seconds to wait between polling batch')
    cl.add_argument('--loglevel', type=int, default=loglevel_default,
                    help='detail for reporting REST queries')
    # poll_count should be zero so each job is polled only once per batch processing loop
    # if poll_count = 1, poll time is not used
    # cl.add_argument('poll_time', type=int, default=20,_default, help='Seconds to wait between polling a single job')
    # cl.add_argument('poll_count', type=int, default=1,_default, help='Number of times to poll each job per batch loop')
    cl.add_argument('fasta_in', type=argparse.FileType('r'))

    return cl.parse_args()


def batch_process(jobs_pending, n_sequence, batch_limit, batch_wait):
    """---------------------------------------------------------------------------------------------
    Poll the submitted jobs until all have finished and return the output. Outputlist is empty if
    1) the batch limit has not been reached
    2) the output is empty (no hits)

    :param jobs_pending: list of Interpro objects
    :param n_sequence: int, total number of sequences submitted
    :param batch_limit: int, number of sequences per batch
    :return: list of string, output
    ---------------------------------------------------------------------------------------------"""
    if (n_sequence % batch_limit):
        return []

    # batch processing loop - wait for all jobs to finish
    # jobs_pending is False when the maximum number of sequences have been submitted
    outputlist = []
    while jobs_pending:
        # poll all pending jobs (ips.status includes poll_time seconds pause). Copy jobs that
        # have not terminated into unfinished, and make this the new list of pending jobs for
        # the next round
        unfinished = []
        time.sleep(batch_wait)

        for i in range(len(jobs_pending)):
            ips = jobs_pending[i]
            if ips.status():
                ips.result()
                if ips.content:
                    # content may be empty if there are no hits
                    outputlist += ips.content.split('\n')

            else:
                unfinished.append(jobs_pending[i])

        jobs_pending = unfinished.copy()

    return outputlist

def poll_all( jobs_pending, batch_wait=61, max_tries=50):
    """---------------------------------------------------------------------------------------------
    Poll the jobs in the jobs_pending list until all have finished. Finished can be
        1) reached maximum number of polling attempts
        2) returned a status other than success or waiting
        3) success

    :param jobs_pending: list of interproscan objects that have been submitted
    :param max_tries: int, maximum number of times to poll
    :return:
    ---------------------------------------------------------------------------------------------"""

    # jobs_pending is False when the maximum number of sequences have been submitted
    outputlist = []
    while jobs_pending:
        # poll all pending jobs (ips.status includes poll_time seconds pause). Copy jobs that
        # have not terminated into unfinished, and make this the new list of pending jobs for
        # the next round
        unfinished = []
        time.sleep(batch_wait)

        for ips in jobs_pending:
            if ips.status() == 'FINISHED':
                ips.result()
                if ips.content:
                    # content may be empty if there are no hits
                    outputlist += ips.content.split('\n')

            else:
                unfinished.append(ips)

        jobs_pending = unfinished.copy()

    return outputlist


# ==================================================================================================
# Main
# ==================================================================================================

args = arguments_get()
args.logfile.write('\ninterpro_batch - interproscan of ORF sequences\n')
args.logfile.write('\tinput ORF file: {}\n'.format(args.fasta_in.name))
args.logfile.write('\tminimum ORF length: {}\n\n'.format(args.minlen))

fasta = Fasta(fh=args.fasta_in)

batch_limit = 2
batch_wait = 30
total_seq = 3

job_list = {}
jobs_pending = []
outputlist = []
n_sequence = 0

while fasta.next():
    if n_sequence >= total_seq:
        # run up to total_seq queries
        break

    if fasta.length() < args.minlen:
        # skip short sequences
        continue

    proseq = fasta.translate()
    proseq.seq = proseq.seq.rstrip('*')
    n_sequence += 1

    # submit jobs
    ips = Interpro(loglevel=1, poll_time=batch_wait, poll_max=1)
    ips.log_fh = args.logfile
    ips.email = 'gribskov@purdue.edu'
    ips.title = 'ORF{}'.format(n_sequence)
    ips.sequence = proseq.format(linelen=60)
    ips.application_select(['Pfam'])
    if ips.run():
        jobs_pending.append(ips)
    else:
        sys.stderr.write('Error - sequence={}\n'.format(fasta.id),flush=True)
        continue

    if n_sequence % batch_limit and n_sequence < total_seq:
        # keep submitting until batch_limit is reached
        continue

    # replace with poll_all
    # outputlist = batch_process(jobs_pending, n_sequence, batch_limit, batch_wait)
    poll_all(jobs_pending)
    for job in outputlist:
        # sys.stdout.write('{}\n'.format(job.content))
        jobs_pending = []

    # if n_sequence > 50:
    #     break

# end of loop over all sequences

outputlist = batch_process(jobs_pending, n_sequence, 1, batch_wait)
for line in outputlist:
    sys.stdout.write('{}\n'.format(line))
    jobs_pending = []

args.logfile.close()
exit(0)

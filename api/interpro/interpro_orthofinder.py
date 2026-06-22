"""=================================================================================================
Run interproscan using the EBI job dispatcher using orthofinder OGs as input. OGs are in FastA
format with all sequences in the OG in a single file

Usage
    intropro_orthofinder.py <OG_directory>

26 December 2018    Michael Gribskov
================================================================================================="""
import sys
import time
from pathlib import Path
import argparse
from sequence.fasta import Fasta
from api.interpro.interpro import Interpro
import textwrap as _textwrap


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
    '''=============================================================================================
    Custom formatter for command line argument help
    ============================================================================================='''

    def _split_lines(self, text, width=60):
        '''-----------------------------------------------------------------------------------------
        Gracefully split lines in command line help. lines are split at 60 characters by default

        :param text: str, text to split
        :param width: int, width to split at
        :return:
        -----------------------------------------------------------------------------------------'''
        text = self._whitespace_matcher.sub(' ', text).strip()
        self._max_help_position = 30
        return _textwrap.wrap(text, width)


# end of class CustomFormatter


def arguments_get():
    '''---------------------------------------------------------------------------------------------
    Set up command line arguments, read from command line, and store in argparse.ArgumentParser
    object, cl

    :return: argparse.ArgumentParser object
    ---------------------------------------------------------------------------------------------'''
    minlen_default = 50
    log_level_default = 1
    batch_limit_default = 30
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
    cl.add_argument('--log_level', type=int, default=log_level_default,
                    help='detail for reporting REST queries')
    cl.add_argument('--og_dir', type=str)

    return cl.parse_args()  # parse_args  reads the command line


def poll_all(joblist, poll_time=61, poll_max=50):
    """---------------------------------------------------------------------------------------------
    Poll the jobs in the jobs_pending list until all have finished. Finished can be
        1) reached maximum number of polling attempts
        2) returned a status other than success or waiting
        3) success

    :param joblist: list of interproscan objects that have been submitted
    :param poll_time: int, seconds to wait between polling
    :param poll_max: int, maximum number of times to poll
    :return: int, number of jobs in list
    ---------------------------------------------------------------------------------------------"""

    not_all_finished = True
    n = 0
    while not_all_finished:
        n += 1
        not_all_finished = False
        time.sleep(poll_time)

        for ips in joblist:
            if ips.status() == 'finished':
                joblist[ips] = 'finished'

            else:
                not_all_finished = True

    return n


def reformat(job):
    """---------------------------------------------------------------------------------------------
    Return a text string with the result of an interproscan job processed with ips.parse_json()
    An example of a callback function for parsing output

    :param job: interpro object, should be a finished job
    :return: string
    ---------------------------------------------------------------------------------------------"""
    str = ''

    parsed = job.parse_json()
    motifs = parsed['motifs']
    go = parsed['go']
    path = parsed['pathway']

    for m in motifs:
        str += '{}\t{}\t{}\n'.format(m['ipr_accession'],
                                     m['src_accession'],
                                     m['description'])
    for g in go:
        str += '{}\t{}\t{}\t{}\n'.format(g, go[g]['name'], go[g]['category'], go[g]['source'])

    for p in path:
        str += '{}\t{}\t{}\n'.format(p, path[p]['name'], path[p]['source'])

    return str


def save_finished(joblist, reformat=None, fh=None, remove=True):
    """---------------------------------------------------------------------------------------------
    Return the output of all finished jobs as a string.
    Reformat is a callback function used to reformat the output.  for instance, interpro.parse_json
    If fh is True, output is written to the filehandle after reformatting.
    If remove is true, jobs are delete from the list after saving

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


def fasta_read(fh, doc=''):
    """---------------------------------------------------------------------------------------------
    Read in the FastA sequences for the orthogroup as a dict of strings. String is the complete
    FastA sequence, key is the ID

    :param fh: filehandle   open for reading
    :param doc: string      extra string to add to docline
    :return: dict           FastA sequence strings
    ---------------------------------------------------------------------------------------------"""
    og_seqs = {}
    for line in fh:
        if line.startswith('>'):
            # for orthofinder input, we expect just a sequence id
            sid = line.rstrip()[1:]
            og_seqs[sid] = f'>{sid} {doc}\n'
        else:
            og_seqs[sid] += line.rstrip('\n*')

    return og_seqs


def fasta_to_batch(fasta, seq_per_batch):
    """---------------------------------------------------------------------------------------------
     make a list of  multi-fasta sequences so that each contains no more than seq_per_batch
     sequences and each partition is as similar as possible in size

    :param fasta: dict              each item is a fasta sequence
    :param seq_per_batch: integer   maximum number of sequences per batch
    :return: list of strings        the sequences to be submitted in each batch
    ---------------------------------------------------------------------------------------------"""
    batchsize = (len(fasta) - 1) // seq_per_batch + 1
    print(f'batchsize: {batchsize}')
    set_n = (len(fasta) + 1) // 2
    nseq = 0
    seq = []
    for id in fasta:
        s = fasta[id]
        if nseq % set_n == 0:
            seq.append('')

        seq[-1] += s + '\n'
        nseq += 1

    return seq


# ==================================================================================================
# Main
# ==================================================================================================

args = arguments_get()
args.logfile.write('\ninterpro_orthofinder - interproscan of Orthofinder OGs\n')
args.logfile.write('\tOG directory: {}\n'.format(args.og_dir))
args.logfile.write('\tminimum ORF length: {}\n\n'.format(args.minlen))

# fasta = Fasta(fh=args.fasta_in)

# The job list keeps track of the ips object that have been created and their current status
# the joblist is a dictionary where the ips object is the key and the value is a status string
joblist = {}

# create a template for the jobs.  The template is an interpro object with the metadata added
template = Interpro(log_level=1)
template.log_fh = args.logfile
template.email = 'gribskov@purdue.edu'
template.application_select(['Pfam', 'Panther', 'SignalP-Euk'])
template.output_select = 'gff'
template.poll_time = 10
template.poll_max = 25
sequence_per_query = 30
simultaneous_jobs = 1

ogfiles = [f.name for f in Path(args.og_dir).iterdir() if f.is_file()]
for f in ogfiles:
    og_f = f'{args.og_dir}/{f}'
    print(f'og file: {og_f}')
    og = open(og_f, 'r')

    # read and store fasta
    fasta = fasta_read(og, f'OG: {f}')
    og.close()

    batch_seq = fasta_to_batch(fasta, sequence_per_query)

    for seq in batch_seq:
        # copy the template and add the sequence information
        ips = template.clone()
        ips.sequence = seq
        ips.jobname = "{f.replace('.fa', '')}"
        ips.title = "{f.replace('.fa', '')}"
        ips.submit()
        joblist[ips] = 'submitted'

        # polling - wait for the job to finish before submitting another
        poll_all(joblist, template.poll_time, template.poll_max)
        save_finished(joblist, reformat, sys.stdout, True)
        exit(1)

    # end of loop over sequence batches

# end of loop over OG groups

args.logfile.close()
exit(0)

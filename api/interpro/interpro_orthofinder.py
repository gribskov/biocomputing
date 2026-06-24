"""=================================================================================================
Run interproscan using the EBI job dispatcher using orthofinder OGs as input. OGs are in FastA
format with all sequences in the OG in a single file

Usage
    intropro_orthofinder.py <OG_directory>

jobmanager_api:JobManager defines the general interface and submitting/polling/retrieving jobs
interpro:InterproscanAPI provides the required functions defined as abstract methods in
    JobManager
interpro:InterproscanQuery holds the information required for queries and methods needed for
    processing results

22 June 2026    Michael Gribskov
================================================================================================="""
import sys
import time
from pathlib import Path
import argparse
from sequence.fasta import Fasta
from api.interpro.interpro import InterproscanQuery, InterproscanAPI
import textwrap as _textwrap


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """=============================================================================================
    Custom formatter for command line argument help
    ============================================================================================="""

    def _split_lines(self, text, width=60):
        """-----------------------------------------------------------------------------------------
        Gracefully split lines in command line help. lines are split at 60 characters by default

        :param text: str, text to split
        :param width: int, width to split at
        :return:
        -----------------------------------------------------------------------------------------"""
        text = self._whitespace_matcher.sub(' ', text).strip()
        self._max_help_position = 30
        return _textwrap.wrap(text, width)


# end of class CustomFormatter


def arguments_get():
    """---------------------------------------------------------------------------------------------
    Set up command line arguments, read from command line, and store in argparse.ArgumentParser
    object, cl

    :return: argparse.ArgumentParser object
    ---------------------------------------------------------------------------------------------"""
    minlen_default = 50  # not implemented
    log_level_default = 1
    sequences_per_query_default = 30
    poll_delay = 60
    simultaneous_jobs_default = 1
    ogdir_default = './data'
    cl = argparse.ArgumentParser(description='Interproscan of ORF sequences',
                                 formatter_class=CustomFormatter)
    cl.add_argument('--logfile', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output file for log information')
    cl.add_argument('-m', '--minlen', type=int, default=minlen_default,
                    help='Minimum length ORF to run')
    cl.add_argument('--sequences_per_query', type=int, default=sequences_per_query_default,
                    help='Number of sequences to submit per batch')
    cl.add_argument('--simultaneous_jobs', type=int, default=simultaneous_jobs_default,
                    help='Number of sequences to submit per batch')
    cl.add_argument('--poll_delay', type=int, default=poll_delay,
                    help='Seconds to wait between polling jobs')
    cl.add_argument('--log_level', type=int, default=log_level_default,
                    help='detail for reporting REST queries')
    cl.add_argument('--ogdir', type=str, default=ogdir_default,
                    help='directory to store orthogroup FastA files')

    return cl.parse_args()  # parse_args  reads the command line


def reformat(job):
    """---------------------------------------------------------------------------------------------
    Return a text string with the result of an interproscan job processed with ips.parse_json()
    An example of a callback function for parsing output

    :param job: interpro object, should be a finished job
    :return: string
    ---------------------------------------------------------------------------------------------"""
    outstr = ''

    parsed = job.parse_json()
    motifs = parsed['motifs']
    go = parsed['go']
    path = parsed['pathway']

    for m in motifs:
        outstr += '{}\t{}\t{}\n'.format(m['ipr_accession'],
                                        m['src_accession'],
                                        m['description'])
    for g in go:
        outstr += '{}\t{}\t{}\t{}\n'.format(g, go[g]['name'], go[g]['category'], go[g]['source'])

    for p in path:
        outstr += '{}\t{}\t{}\n'.format(p, path[p]['name'], path[p]['source'])

    return outstr


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
    text = None
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
    sid = ''
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
    for sid in fasta:
        s = fasta[sid]
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
args.logfile.write('\tOG directory: {}\n'.format(args.ogdir))
args.logfile.write('\tminimum ORF length: {}\n\n'.format(args.minlen))


# manager handles the specifics of submitting jobs, polling, and retrieving results. Manager
# is reused for each query (which is an InterproscanQuery object)
manager = InterproscanAPI()

# query values that do not change for each query
# TODO validate applications and output
constants = {'url': u'https://www.ebi.ac.uk/Tools/services/rest/iprscan6/',
             'program': 'iprscan6',
             'email': 'gribskov@purdue.edu',
             'appl': ['SUPERFAMILY', 'Pfam', 'Panther', 'SignalP-Euk'],
             'resultType': 'gff3',
             'goterms': True,
             'pathways': False,
             'sequence': '',
             'stype': 'p',
             'title': '' }

ogfiles = [f.name for f in Path(args.ogdir).iterdir() if f.is_file()]
for f in ogfiles:
    og_f = f'{args.ogdir}/{f}'
    print(f'og file: {og_f}')
    og = open(og_f, 'r')

    # read and store fasta
    fasta = fasta_read(og, f'OG: {f}')
    og.close()

    batch_seq = fasta_to_batch(fasta, args.sequences_per_query)
    batches = len(batch_seq)

    batch_num = 0
    for seq in batch_seq:
        batch_num += 1
        query = InterproscanQuery(constants)
        query.parameters['sequence'] = seq
        query.parameters['jobname'] = f"{f.replace('.fa', '')}"
        if batches > 1:
            query.parameters['jobname'] += f'_{batch_num}'
        query.parameters['title'] = f"{f.replace('.fa', '')}"
        query.validate(['appl', 'resultType'])

        # Iprscan service says to wait for the job to finish before submitting another
        print(f'\nsubmitting {query.parameters["jobname"]}')
        manager.submit(query, show_query=True)
        manager.joblist.append(query)
        manager.poll_all()
        print(f'finished {query.parameters["jobname"]}')
        manager.result_all()
        manager.save_all()
        print(f'writing {query.parameters["jobname"]}')
        manager.joblist.remove(query)
        exit(1)

    # end of loop over sequence batches

# end of loop over OG groups

args.logfile.close()
exit(0)

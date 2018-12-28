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


def arguments_get():
    minlen_default = 80
    cl = argparse.ArgumentParser(description='Interproscan of ORF sequences')
    cl.add_argument('-m', '--minlen', type=int, default=minlen_default,
                    help='minimum length ORF to run')
    cl.add_argument('fasta_in', type=argparse.FileType('r'))

    return cl.parse_args()


# ==================================================================================================
# Main
# # ==================================================================================================
# batch_limit = 20
# batch_wait = 300
batch_limit = 4
batch_wait = 30

args = arguments_get()
sys.stderr.write('\ninterpro_batch - interproscan of ORF sequences\n')
sys.stderr.write('\tinput ORF file: {}\n'.format(args.fasta_in.name))
sys.stderr.write('\tminimum ORF length: {}\n\n'.format(args.minlen))

fasta = Fasta(fh=args.fasta_in)
n_sequence = 0

job_list = {}
jobs_pending = []
while fasta.next():
    if fasta.length() >= args.minlen:
        n_sequence += 1
        ips = Interpro(loglevel=1, poll_time=20, poll_count=1)
        ips.email = 'gribskov@purdue.edu'
        ips.title = 'ORF{}'.format(n_sequence)
        ips.sequence = fasta.format(linelen=60)
        if not ips.run():
            sys.stderr.write('Error - sequence={}\n'.format(fasta.id))
            continue
        jobs_pending.append(ips)
    else:
        continue

    if not (n_sequence % batch_limit):
        # False when the maximum number of sequences have been submitted
        # wait for all sequences to finish
        while jobs_pending:
            unfinished = []
            for i in range(len(jobs_pending)):
                ips = jobs_pending[i]
                if ips.status():
                    ips.result()
                    sys.stdout.write(ips.content)

                else:
                    unfinished.append(jobs_pending[i])

            jobs_pending = unfinished.copy()
            time.sleep(batch_wait)

    if n_sequence > 3:
        break

# end of loop over all sequences

exit(0)
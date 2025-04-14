"""-------------------------------------------------------------------------------------------------
select sequences from a Fasta file using
    list of IDs
    minimum length

usage:
    # select sequences from a single fasta file based on a list of desired sequences and write to
    # stdout, format fasta output with 60 character lines
    fasta_select in.fa --line 60 --list select.txt

    # select sequences from a single fasta file based minimum length and write to
    # stdout, format fasta output with 100 character lines
    fasta_select in.fa --line 100 --min 600

    # select sequences from a single fasta file skipping blast matches, with min length cutoff
    fasta_select in.fa --blast search.blastn --evalue 1e-20 --min 2000

    # select sequences from a multiple fasta files based minimum length and write to
    # stdout, format fasta output with 100 character lines
    fasta_select genome*.fa --line 100 --min 600

    # select sequences from a multiple fasta files based minimum length and write to
    # individual files with a suffix. Output files are in the directory where the program is run
    genome12.fa -> genome12.fa.select
    fasta_select genome*.fa --line 100 --min 600 --out .select

    # select sequences from a multiple fasta files based minimum length and write to
    # individual files with a suffix. Trim the documentation line to remove some information (
    # such as the trinity path)
    fasta_select genome*.fa --line 100 --min 600 --out .select --trim "path=.*"

2 September 2018    Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import glob
import os.path
import sys
import time
import re
import argparse
from sequence.fasta import Fasta


def setup_argparse():
    """---------------------------------------------------------------------------------------------
    Setup parsing of command line arguments with argparse

    :return: argparse ArgumentParser
    ---------------------------------------------------------------------------------------------"""
    # defaults
    # out_default = sys.stdout
    minlen_default = 0
    linelen_default = 100
    blast_default = ''
    list_default = 'list.default'
    evalue_default = 1e-5

    commandline = argparse.ArgumentParser(
        description='Select sequences from a multiple-sequence FastA file'
        )

    # the input file is a string because it may be a wildcard
    commandline.add_argument('input_filename',
                             help='FastA file to select from. Wildcards are allowed.',
                             type=str)

    commandline.add_argument('--outsuffix',
                             help='Suffix for output file. Output will be STDOUT if not provided.',
                             type=str,
                             default='',
                             )

    commandline.add_argument('--list',
                             help='list of sequence names to select (all).',
                             type=str,
                             default=list_default
                             )

    commandline.add_argument('--trim',
                             help=f'Regular expression for trimming doc line(< none >).',
                             default='')

    commandline.add_argument('--minlen',
                             help=f'Minimum length for selected sequences ({minlen_default}).',
                             type=int,
                             default=minlen_default)

    commandline.add_argument('--linelen',
                             help=f'Line length for FastA format ({linelen_default}).',
                             type=int,
                             default=linelen_default)

    commandline.add_argument('--blast',
                             help=f'Filter based on blast search ({blast_default}).',
                             type=str,
                             default=blast_default)

    commandline.add_argument('--evalue',
                             help=f'Maximum E-value for blast filtering ({evalue_default}).',
                             type=float,
                             default=evalue_default)

    return report_args(commandline.parse_args())


def report_args(args):
    """---------------------------------------------------------------------------------------------
    Report command line parameters and other useful information to STDERR

    :param args: dict, arguments from setup_argparse()
    :return: args
    ---------------------------------------------------------------------------------------------"""
    now = time.asctime(time.localtime(time.time()))
    sys.stderr.write(f'fasta.select {now}\n')
    sys.stderr.write(f'\tinput: {args.input_filename}\n')
    if args.outsuffix:
        sys.stderr.write(f'\toutput: <input_file>.{args.outsuffix}\n')
    else:
        sys.stderr.write(f'\toutput file: STDOUT\n')
    sys.stderr.write(f'\tFastA line length: {args.linelen}\n')
    if args.minlen:
        sys.stderr.write(f'\tMinimum sequence length: {args.minlen}\n')
    if args.list:
        sys.stderr.write(f'\tSequence ID list: {args.list}\n')
    if args.trim:
        sys.stderr.write(f'\tDocumentation trimming regex: {args.trim}\n')
    if args.blast:
        sys.stderr.write(f'\tFilter using Blast file: {args.blast}\n')
        sys.stderr.write(f'\tE-value cutoff: : {args.evalue}\n')
    sys.stderr.write('\n')

    return args


def opensafe(filename, mode, die_on_error=False):
    """---------------------------------------------------------------------------------------------
    open a file with error checking.  if die_on_error is true exit with status=1 on failure.

    :param filename: str
    :param mode: str
    :param die_on_error: boolean
    :return: file handle
    ---------------------------------------------------------------------------------------------"""
    file = None
    try:
        file = open(filename, mode)
    except (OSError, IOError):
        sys.stderr.write(f'opensafe: unable to open ({filename}) in mode({mode})')
        if die_on_error:
            exit(1)

    return file


def get_id_list(args):
    """---------------------------------------------------------------------------------------------
    Read a file of IDs and store in a list.  Only the first token is used
    Lines beginning ! or # are skipped, if > is the first character it is removed

    :param args: dict, command line arguments from setup_argparse()
    :return: list of strings, IDs of desired sequences
    ---------------------------------------------------------------------------------------------"""
    if not args.list:
        # blank list if no idlist is given
        return []

    # open and read list
    idfile = opensafe(args.list, 'r', die_on_error=True)

    nid = 0
    idlist = []
    for line in idfile:
        if line.startswith(('#', '!')):
            continue
        line = line.lstrip('>')
        if line.rstrip():
            # remove double quotes and split
            id = line.replace('"', '').split()[0]
            if id == 'baseMean':
                # this looks like it deals with input from DeSEQ2
                continue
            if id not in idlist:
                # make sure IDs are unique in list
                idlist.append(id)
                nid += 1

    sys.stderr.write(f'{nid} sequence IDs read from {args.list}\n')
    return idlist


def blast_filter(blastfilename, evalue):
    """---------------------------------------------------------------------------------------------
    Based on a blast search, make a list of query sequences that have matches <= evalue. Use this to
    remove contaminants such as matches to host sequence in a pathogen dataset

    :param blastfilename: string        path to openable blast result (format=7)
    :param evalue: float                maximum evalue
    :return: list                       query sequence IDs
    ---------------------------------------------------------------------------------------------"""
    blast = opensafe(blastfilename, 'r')

    seq_n = 0
    select_n = 0
    select = []
    for line in blast:
        if line.startswith('# Query:'):
            seq_n += 1
            continue
        elif line.startswith('#'):
            continue
        else:
            field = line.split()
            if field[0] in select:
                continue
            if float(field[10]) > evalue:
                continue
            select.append(field[0])
            select_n += 1

        if seq_n > 100:
            break

    blast.close()
    return select


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    args = setup_argparse()

    # set up regular expression for trimming documentation
    trim = ''
    if args.trim:
        trim = re.compile(args.trim)

    # idlist, idlist will be an emtpy list if none is provided
    idlist = get_id_list(args)

    blastfilter = False
    blastlist = []
    if args.blast:
        blastfilter = True
        evalue = args.evalue
        blastlist = blast_filter(args.blast, args.evalue)
        sys.stderr.write(f'sequences filtered based on {args.blast}: {len(blastlist)}\n')

    # read the sequences and store all that match the IDs
    # duplicates in sequence files will be stored twice
    n_match = {}  # per file number of sequences in list
    n_notmatch = {}  # per file number of sequences not in list
    n_sequence = {}  # per file number of sequences
    n_found = {}  # per ID, number of times found in all files
    n_file = 0
    n_total = 0
    n_written = 0
    out = sys.stdout
    for fastafile in glob.glob(args.input_filename):
        fasta = Fasta()
        fasta.open(fastafile)
        if args.outsuffix:
            outfile = os.path.basename(fastafile) + f'{args.outsuffix}'
            out = opensafe(outfile, 'w')
            if not out:
                # if file can't be opened use stdout
                out = sys.stdout

        n_sequence[fastafile] = 0
        n_match[fastafile] = 0
        n_notmatch[fastafile] = 0
        n_file += 1

        while fasta.next():
            n_sequence[fastafile] += 1
            n_total += 1
            if n_total % 1000 == 0:
                sys.stderr.write('.')
            if n_total % 50000 == 0:
                sys.stderr.write(f'{n_total}\n')

            # if blastfilter is true, check if this sequence should be skipped
            if blastfilter and fasta.id in blastlist:
                continue

            if fasta.id in idlist or not idlist:
                # desired selected sequences
                if args.trim:
                    fasta.trimDocByRegex(trim)
                seqlen = len(fasta.seq)
                if args.minlen and seqlen < args.minlen:
                    # skip sequences shorter than minimum length, if specified
                    continue

                out.write('{}\n'.format(fasta.format(linelen=args.linelen)))
                n_written += 1
                n_match[fastafile] += 1
                if fasta.id in n_found:
                    n_found[fasta.id] += 1
                else:
                    n_found[fasta.id] = 1

            else:
                # not selected sequence
                n_notmatch[fastafile] += 1

    sys.stderr.write('\nfiles read: {}\n'.format(n_file))
    sys.stderr.write('total sequences read: {}\n'.format(n_total))
    sys.stderr.write('total sequences written: {}\n'.format(n_written))

    if length(n_sequence) > 1:
        sys.stderr.write('\nPer file\n')
        for fastafile in n_sequence:
            sys.stderr.write('{}\n'.format(fastafile))
            sys.stderr.write('\tsequences: {}\n'.format(n_sequence[fastafile]))
            sys.stderr.write('\tsequences matched: {}\n'.format(n_match[fastafile]))

exit(0)

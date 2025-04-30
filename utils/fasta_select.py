"""-------------------------------------------------------------------------------------------------
select sequences from a Fasta file using
    list of IDs or blast result
    minimum length
    evalue cutoff (blast result only)
    filter input ID list to match IDs in fasta file
    filter fasta documentation (to remove excessively long documentation lines)
    use options --idcol and --evalue to select the columns for the IDs and evalues respectively

usage:
    # select sequences from a single fasta file based on a list of desired sequences and write to
    # stdout, format fasta output with 60 character lines
    python fasta_select.py --line 60 --list select.txt in.fa 

    # select sequences from a single fasta file based minimum length (=600) and write to
    # stdout, format fasta output with 100 character lines
    python fasta_select.py --line 100 --min 600 in.fa 

    # select sequences from a single fasta file based on a blast result. select blast matches with
    # E <= 1e-20 (column 10) and length >= 2000
    python fasta_select.py --format blast --blast search.blastn --evalue 1e-20,10 --min 2000 in.fa

    # select sequences from a single fasta file based on a blast result. select blast matches with
    # E <= 1e-20 and length >= 2000s. Remove the prefix TRINITY_ from the IDs in the blast file
    python fasta_select.py --format blast --blast search.blastn --trimid 'TRINITY_'
                           --evalue 1e-20,10 --min 2000 in.fa

    # select sequences from a multiple fasta files based minimum length and write to stdout, format
    # fasta output with 100 character lines
    python fasta_select.py  --line 100 --min 600 genome*.fa

    # select sequences from a multiple fasta files based minimum length and write to individual
    # files with a suffix (.select), genome12.fa -> genome12.fa.select. Output files are written to
    # the directory where the program is run
    python fasta_select.py --line 100 --min 600 --out .select genome*.fa 

    # select sequences from multiple fasta files based minimum length and write to individual files
    # with the suffix .select. Trim the documentation line to remove some information (such as the
    # trinity assembly path)
    python fasta_select.py --line 100 --min 600 --out .select --trimdoc "path=.*" genome*.fa 

2 September 2018    Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import argparse
import glob
import os.path
import re
import sys
import time

from sequence.fasta import Fasta


def setup_argparse():
    """---------------------------------------------------------------------------------------------
    Setup parsing of command line arguments with argparse

    :return: argparse ArgumentParser
    ---------------------------------------------------------------------------------------------"""
    # defaults
    # out_default = sys.stdout
    idcol_default = 0
    trinity_level_default = 0
    minlen_default = 0
    linelen_default = 100
    list_default = 'list.txt'
    format_default = 'list'
    evalue_col = 10
    evalue_max = 1e-5
    evalue_default = f'{evalue_max},{evalue_col}'

    commandline = argparse.ArgumentParser(
        description='Select sequences from a multiple-sequence FastA file'
        )

    # the input file is a string because it may be a wildcard
    commandline.add_argument('input_filename',
                             help='FastA file to select from. Wildcards are allowed.',
                             type=str)

    commandline.add_argument('--outsuffix',
                             help='Suffix for output file (same prefix as list input). Output will '
                                  'be STDOUT if not provided.',
                             type=str,
                             default='',
                             )

    commandline.add_argument('--format',
                             help='Format for sequence list (blast|list)',
                             type=str,
                             default=format_default
                             )

    commandline.add_argument('--list',
                             help='list of sequence names to select.',
                             type=str,
                             default=list_default
                             )

    commandline.add_argument('--idcol',
                             help=f'Column of sequence ID in selection list ({idcol_default})',
                             type=int,
                             default=idcol_default
                             )

    commandline.add_argument('--trinity',
                             help=f'Trim database ids at specified level, zero means do not trim'
                                  f'({trinity_level_default}).',
                             default='trinity_level_default')

    commandline.add_argument('--evalue',
                             help=f'Filter blast search by e-value: max_evalue,column ({evalue_default}).',
                             type=str,
                             default=evalue_default)

    commandline.add_argument('--trimid',
                             help=f'Regular expression for trimming IDs(< none >).',
                             default='')

    commandline.add_argument('--trimdoc',
                             help=f'Regular expression for trimming doc line(< none >).',
                             default='')

    commandline.add_argument('--minlen',
                             help=f'Minimum length for selected sequences ({minlen_default}).',
                             type=int,
                             default=minlen_default)

    commandline.add_argument('--linelen',
                             help=f'Line length for FastA format output ({linelen_default}).',
                             type=int,
                             default=linelen_default)

    cl = commandline.parse_args()

    # check that list format is allowed
    if cl.format == 'list':
        cl.listread = list_read
    elif cl.format == 'blast':
        cl.listread = blast_read
    else:
        sys.stderr.write('Unknown input list format ({cl.format}. Select from list|blast')
        exit(1)

    # get e-value cutoff and column, only meaningful for blast search as selection input
    field = cl.evalue.split(',')
    cl.evalue = float(field[0])
    if len(field) > 1:
        cl.evalue_col = int(field[1])

    # setup regular expressions for modifying selected IDs and fasta documentation
    if cl.trimid:
        cl.trimid = re.compile(cl.trimid.replace("'", ""))
    if cl.trimdoc:
        cl.trimdoc = re.compile(cl.trimdoc.replace("'", ""))

    return report_args(cl)


def report_args(args):
    """---------------------------------------------------------------------------------------------
    Report command line parameters and other useful information to STDERR

    :param args: dict, arguments from setup_argparse()
    :return: args
    ---------------------------------------------------------------------------------------------"""
    now = time.asctime(time.localtime(time.time()))
    sys.stderr.write(f'fasta_select.py {now}\n')
    sys.stderr.write(f'\tFastA filet: {args.input_filename}\n')
    if args.format == 'list':
        sys.stderr.write(f'\tSelection IDs read from list: {args.list}\n')
    if args.format == 'blast':
        sys.stderr.write(f'\tSelection IDs read from Blast file: {args.list}\n')
        sys.stderr.write(f'\tE-value cutoff: {args.evalue}\n')

    if args.outsuffix:
        sys.stderr.write(f'\tOutput: <input_file>.{args.outsuffix}\n')
    else:
        sys.stderr.write(f'\tOutput file: STDOUT\n')

    sys.stderr.write(f'\tFastA line length: {args.linelen}\n')
    if args.minlen:
        sys.stderr.write(f'\tMinimum sequence length: {args.minlen}\n')

    if args.trimdoc:
        sys.stderr.write(f'\tDocumentation trimming regex: {args.trimdoc}\n')
    if args.trimid:
        sys.stderr.write(f'\tID trimming regex: {args.trimid}\n')

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


def list_read(args, verbose=True):
    """---------------------------------------------------------------------------------------------
    Read a file of IDs and store in a list.  Only the first token is used
    Lines beginning ! or # are skipped, if > is the first character it is removed

    :param args: dict           command line arguments from setup_argparse()
    :param verbose: boolean     write a report to stderr
    :return: list               of strings, IDs of desired sequences
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

        line = line.lstrip('>').rstrip()
        # remove double quotes and split
        sid = line.replace('"', '').split()[0]
        if sid == 'baseMean':
            # this looks like it deals with input from DeSEQ2
            continue

        if sid not in idlist:
            # make sure IDs are unique in list
            idlist.append(sid)
            nid += 1

    if verbose:
        sys.stderr.write(f'{nid} sequence IDs read from {args.list}\n')
    idfile.close()
    return idlist


def idlist_filter(args, idlist):
    """---------------------------------------------------------------------------------------------
    Filter the IDs in idlist using the regular expression in

    :param args: namespace      command line options from setup_argparse()
    :param idlist: list         IDs from input selection list
    :return: None
    ---------------------------------------------------------------------------------------------"""
    if not args.trimid:
        # no trimming regex
        return

    trimre = args.trimid
    for i in range(len(idlist)):
        idlist[i] = trimre.sub('', idlist[i])

    return None

def trinity_filter(idstr, level):
    """---------------------------------------------------------------------------------------------
    return a trinity id truncated at level. if level==0 return the original ID
    level   id
        4   TRINITY_DN11933_c17_g1_i1
        3   TRINITY_DN11933_c17_g1
        2   TRINITY_DN11933_c17
        1   TRINITY_DN11933

    :param idstr: str   trinity ID
    :param level: int   truncation level
    :return: str        truncated ID
    ---------------------------------------------------------------------------------------------"""
    field = idstr.split('_')
    return '_'.join(field[:level+1])


def blast_read(args, verbose=True):
    """---------------------------------------------------------------------------------------------
    Based on a blast search, make a list of query sequences that have matches <= evalue. Use this to
    remove contaminants such as matches to host sequence in a pathogen dataset
    file to be read is given on command line and saved in args.list

    :param args: namespace              argument list from command line
    :param verbose: boolean            write a report to stderr
    :return: list                       query sequence IDs
    ---------------------------------------------------------------------------------------------"""
    blast = opensafe(args.list, 'r')
    idcol = args.idcol
    evalue = args.evalue
    evaluecol = args.evalue_col

    line_n = 0
    seq_n = 0
    select_n = 0
    select = []
    sid_old = ''
    for line in blast:
        if line.startswith('#'):
            continue
        else:
            line_n += 1
            field = line.split()
            sid = field[idcol]
            if sid != sid_old:
                seq_n += 1
                sid_old = sid

            if sid in select:
                # should not be necessary, but maybe there are several results being combined
                continue

            # new id
            if float(field[evaluecol]) > evalue:
                continue

            select.append(field[idcol])
            select_n += 1

    if verbose:
        sys.stderr.write(f'blast records processed: {line_n} \n')
        sys.stderr.write(f'sequence queries found: {seq_n} \n')

    blast.close()
    return select


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    args = setup_argparse()

    # read list of IDs to extract from fasta file(s)
    idlist = args.listread(args)
    idlist_filter(args, idlist)
    sys.stderr.write(f'{len(idlist)} selection IDs read from {args.list}\n')

    # read the sequences and store all that match the IDs in idlist
    # duplicates in sequence files will be stored twice, this allows sequences that differ in
    # different fasta files to be retained
    n_match = {}  # per file number of sequences in list
    n_notmatch = {}  # per file number of sequences not in list
    n_sequence = {}  # per file number of sequences examined
    n_found = {}  # per ID, number of times found in all files
    n_file = 0
    n_total = 0
    n_written = 0
    out = sys.stdout
    for fastafile in glob.glob(args.input_filename):
        # outer loop selects all fasta files that match the wildcard input
        fasta = Fasta()
        fasta.open(fastafile)
        if args.outsuffix:
            # create output file with the same name as the input with a new suffix
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
            # if blastfilter and fasta.id in blastlist:
            #     continue

            id = trinity_filter(fasta.id, args.trinity_level)
            if id in idlist or not idlist:
                # desired selected sequences
                if args.trimdoc:
                    fasta.trimDocByRegex(args.trimdoc)
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

    if len(n_sequence) > 1:
        sys.stderr.write('\nPer file\n')
        for fastafile in n_sequence:
            sys.stderr.write('{}\n'.format(fastafile))
            sys.stderr.write('\tsequences: {}\n'.format(n_sequence[fastafile]))
            sys.stderr.write('\tsequences matched: {}\n'.format(n_match[fastafile]))

exit(0)

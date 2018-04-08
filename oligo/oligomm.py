"""=================================================================================================
Build oligos that are unlikely in a sequence based on kmer frequencies. Infrequent kmers are
overlapped to probabilistically generate longer kmers unlikely to occur in the sequence.

TODO: check for duplicates?
TODO: bias towards more extreme AT/GC content?
================================================================================================="""
import sys
import argparse
from math import log10 as log

from kmer import Kmer


def commandLine(default):
    """---------------------------------------------------------------------------------------------
    Get command line options

    :return: ArgumentParser.parse_args()
    ---------------------------------------------------------------------------------------------"""
    commandline = argparse.ArgumentParser(
        description='calculate infrequent oligos based on kmer frequencies in an input'
                    'Fasta file. Oligos are written Fasta format to stdout.'

    )
    commandline.add_argument('--fasta',
                             help='FastA file to analyze',
                             type=argparse.FileType('r'),
                             )

    commandline.add_argument('--kmer',
                             help='size of kmer words',
                             type=int,
                             default=str(default['kmer'])
                             )

    commandline.add_argument('--overlap',
                             help='overlap between kmer words',
                             type=int,
                             default=str(default['overlap'])
                             )

    commandline.add_argument('--noligo',
                             help='number of oligos to generate',
                             type=int,
                             default=str(default['noligo'])
                             )

    commandline.add_argument('--loligo',
                             help='length of oligos to generate',
                             type=int,
                             default=str(default['loligo'])
                             )
    commandline.add_argument('--oligo',
                             help='output oligo file',
                             type=argparse.FileType('w')
                             )

    commandline.add_argument('--table',
                             help='precalculated kmer table',
                             type=argparse.FileType('r'),
                             )

    commandline.add_argument('--output',
                             help='output kmer table',
                             type=argparse.FileType('w')
                             )

    return commandline.parse_args()


def printParameters(clargs):
    """---------------------------------------------------------------------------------------------
    Print a report on command line parameters to stderr.

    :param clargs: command line arguments from argparse.parse_args()
    :return: True
    ---------------------------------------------------------------------------------------------"""
    sys.stderr.write('\n')
    if clargs.fasta or clargs.table:
        for k in clargs.__dict__:
            if k.startswith('__') or clargs.__dict__[k] == None:
                continue
            if k in ('fasta', 'table', 'oligo', 'output'):
                sys.stderr.write('{}: {}\n'.format(k, clargs.__dict__[k].name))
            else:
                sys.stderr.write('{}: {}\n'.format(k, clargs.__dict__[k]))

    else:
        sys.stderr.write(
            'Input error: fasta (--fasta) table or kmer table (--table) must be provided\n')
        exit(1)

        sys.stderr.write('\n')
    return True


def fractionGC(seq, alphabet):
    """---------------------------------------------------------------------------------------------
    Calculate fraction GC in a sequence. non ACGT letters are ignored

    :param seq:
    :return: fraction GC
    ---------------------------------------------------------------------------------------------"""
    count = {x: 0 for x in alphabet}
    all = 0
    for base in seq:
        count[base] += 1

        all += 1

    return (count['G'] + count['C']) / all


# ==================================================================================================
# main
# ==================================================================================================

if __name__ == '__main__':

    default = {'kmer': 8, 'noligo': 1000, 'loligo': 50, 'overlap': 3}
    cl = commandLine(default)
    printParameters(cl)

    kmer = Kmer(cl.kmer)

    nwords = 0
    infile = ''
    if cl.table:
        # read kmers from file
        nwords = kmer.tableRead(cl.table)
        sys.stderr.write(
            '\n{} kmers of size {} read from {}\n'.format(kmer.total, kmer.k, cl.table.name))
        infile = cl.table.name
    else:
        # calculate kmers from fasta
        # convert word counts to probabilities and write out
        nwords = kmer.fromFasta(cl.fasta)
        pmin, pmax = kmer.updateProb()
        infile = cl.fasta.name

    # weight counts - ad hoc weighting function
    # w = bias**(logP - logPmax)
    # the larger the bias, the more the weights favor infrequent words, save this for future ref

    kmer.weightExp(-2.0)
    kmer.weightNormalize()
    sys.stderr.write('\n    pmin pmax: {:10.3g}\t{:10.3g}\n'.format(kmer.pmin, kmer.pmax))
    sys.stderr.write('log pmin pmax: {:10.3g}\t{:10.3g}\n'.format(log(kmer.pmin), log(kmer.pmax)))
    sys.stderr.write('    wmin wmax: {:10.3g}\t{:10.3g}\n'.format(kmer.wmin, kmer.wmax))
    sys.stderr.write('\n')

    if cl.output:
        nwritten = kmer.tableWrite(cl.output)
        sys.stderr.write('{} kmers written to {}\n'.format(nwritten, cl.output.name))

    # calculate the transitions between overlapping words: 'list' is the list of overlapping full
    # words, 'wcum' has the threshold value for each word (cumulative weight)

    owords = {}
    for word in kmer.list:

        lap = word[:cl.overlap]
        if lap in owords:
            omer = owords[lap]
        else:
            omer = Kmer(k=cl.overlap, setup=False)
            owords[lap] = omer

        omer.list.append(word)
        omer.add('count', word, kmer.count[word])
        omer.add('p', word, kmer.p[word])
        omer.add('w', word, kmer.w[word])
        omer.wmin = min(omer.wmin, omer.w[word])
        omer.wmax = max(omer.wmax, omer.w[word])
        omer.pmin = min(omer.pmin, omer.p[word])
        omer.pmax = max(omer.pmax, omer.p[word])

    for o in owords:
        owords[o].weightNormalize()

    # construct cl.noligo oligos
    out = sys.stdout
    if cl.oligo is not None:
        out = cl.oligo

    for noligo in range(cl.noligo):
        # select a weighted start point
        oligo, weight = kmer.randomByWeight()
        while len(oligo) < cl.loligo:
            # select an overlapping next word
            lap = oligo[-cl.overlap:]
            next, weight = owords[lap].randomByWeight()
            # sys.stdout.write('    {}  {}  {}\n'.format(lap, next, oligo))
            oligo += next[cl.overlap:]

        gc = fractionGC(oligo, kmer.alphabet)
        out.write('>n{} GC={:.3}  Len={}\n'.format(noligo, gc, len(oligo)))
        out.write('{}\n'.format(oligo))

    sys.stderr.write('{} oligos written to {}\n'.format(cl.noligo, out.name))

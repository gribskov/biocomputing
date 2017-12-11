'''
Trinity
Read multiple trinity transcript file in fasta format
'''
from sequence.fasta import Fasta
import re

lenre = re.compile('len=(\d+)')
pathre = re.compile('.*path=\[([^\]]*)\]')


def getLen(line):
    '''-----------------------------------------------------------------------------------------------------------------
    get the sequence length from the documentation
    :return: length fread from length field in documentation
    -----------------------------------------------------------------------------------------------------------------'''
    return lenre.match(line).group(1)


def getPath(line):
    '''-----------------------------------------------------------------------------------------------------------------
    The path describes how the predicted trasncript is built from segments
    :return: list of path components from documentation
    -----------------------------------------------------------------------------------------------------------------'''
    path = pathre.match(line).group(1)
    plist = path.split(' ')

    return plist


def trinityID(id):
    '''-----------------------------------------------------------------------------------------------------------------
    Breakdown the trinity ID string to give the separate parts of the ID
    Cluster,  component, gene and isoform
    :return: cluster,  component, gene, isoform
     ----------------------------------------------------------------------------------------------------------------'''
    cluster, component, gene, isoform = re.compile(r'>*TRINITY_([^_]+)_c(\d+)_g(\d+)_i(\d+)').match(id).groups()
    return cluster, component, gene, isoform


if __name__ == '__main__':
    trinity = Fasta()

    file = r'C:\Users\gribs\Dropbox\rice\Trinity.fasta'
    trinity.fh = open(file, 'r')

    nseq = 0
    while trinity.next():
        nseq += 1

        cluster, component, gene, isoform = trinityID(trinity.id)
        nn = getPath(trinity.doc)

        if nseq > 10: break
        break

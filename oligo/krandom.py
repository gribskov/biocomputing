"""=================================================================================================
Exact kmer preserving randomization

Michael Gribskov     15 December 2025
================================================================================================="""
import random
from collections import defaultdict


class Word:
    """=============================================================================================

    Synopsis
        (None yet)

    Michael Gribskov
    ============================================================================================="""
    k = 3

    def __init__(self, label=''):
        """-----------------------------------------------------------------------------------------
        word constructor
        -----------------------------------------------------------------------------------------"""
        self.label = ''
        self.out = []
        self.next = None
        self.in_tree = False

        if label:
            self.label = label

def randomtreewithroot(nodelist,rnode):
    """---------------------------------------------------------------------------------------------
    pseudocode in
    Jiang, M., Anderson, J., Gillespie, J. et al. uShuffle: A useful tool for shuffling
    biological sequences while preserving the k-let counts. BMC Bioinformatics 9, 192 (2008).
    https://doi.org/10.1186/1471-2105-9-192

    :return:
    ---------------------------------------------------------------------------------------------"""

    for n in nodelist:
        # mark all nodes as available (not in_tree)
        i = nodelist[n]
        i.in_tree = False

    rnode.next = None
    rnode.in_tree = True

    for n in nodelist:
        u = nodelist[n]
        if u.in_tree:
            continue

        start = u

        while not u.in_tree:
            u.next = random.choice(u.out)
            dumptree(u)
            u = u.next

        u = start
        while not u.in_tree:
            u.in_tree = True
            u = u.next

    return True

def dumptree(node):
    """---------------------------------------------------------------------------------------------
    Trace a tree of nodes by their next pointer starting at node

    :param node: Word       Root node of tree
    :return: True
    ---------------------------------------------------------------------------------------------"""
    while node.next:
        nextnode = node.next
        print(f'node:{node.label} => \t{nextnode.label}',end='\t')
        for i in nextnode.out:
            print(f'{i.label}', end=' ')
        print()
        node = node.next

    return True

# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    k = 3
    alpha = 'ACGT'
    len = 20
    sequence = ''
    for l in range(len):
        base = random.choice(alpha)
        sequence += base
        # print(f'{base}:{sequence}')

    sequence = ('TTGATGGGAACGTATGTGAT')
    print(sequence)

    # Generate De Bruijn graph
    words = defaultdict(Word)
    prev = ''
    for pos in range(0, len - k + 1):
        kmer = sequence[pos:pos + k]
        # print(kmer)

        words[kmer].label = kmer
        # overlap[kmer[1:-1]].out.append(words[kmer])
        if prev:
            words[prev].out.append(words[kmer])
            # overlap[prev[1:-1]].out.append(words[kmer])

        prev = kmer

    # print(words.keys())
    # wordkeys = list(words.keys())
    # randomkey = random.choice(wordkeys)
    # rnode = words[randomkey]
    # end node
    end = words[sequence[-k:]]
    randomtreewithroot(words,end)
    while rnode.next:
        print(rnode.label)
        rnode = rnode.next

    exit(0)

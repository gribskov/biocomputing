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
        self.in_tree = False

        if label:
            self.label = label

def randomtreewithroot(nodelist):
    """-----------------------------------------------------------------------------------------
    pseudocode in
    Jiang, M., Anderson, J., Gillespie, J. et al. uShuffle: A useful tool for shuffling
    biological sequences while preserving the k-let counts. BMC Bioinformatics 9, 192 (2008).
    https://doi.org/10.1186/1471-2105-9-192

    :return:
    -----------------------------------------------------------------------------------------"""

    for i in nodelist:
        i.in_tree = False

    next[r] = None
    in_tree[r] = True
    for i in nodelist:
        u = i

        while not u.in_tree:
            next[u] = randomsuccessor(u)
            u = next[u]


    u = i
    while not u.in_tree:
        u.in_tree = True
        u = next[u]

    return next

def randomsuccessor(self):
    return
        
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
        print(f'{base}:{sequence}')

    # overlap = defaultdict(Word)
    words = defaultdict(Word)
    prev = ''
    for pos in range(0, len - k + 1):
        kmer = sequence[pos:pos + k]
        print(kmer)

        words[kmer].label = kmer
        # overlap[kmer[1:-1]].out.append(words[kmer])
        if prev:
            words[prev].out.append(words[kmer])
            # overlap[prev[1:-1]].out.append(words[kmer])

    exit(0)

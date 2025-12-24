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

# another try
class DBG:
    """

    """

    def __init__(self, sequence='', k=3):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.k = k
        self.sequence = sequence
        self.nodes = set()
        self.root = ''
        self.adj = defaultdict(list)
        self.last_exit = defaultdict(str)

        if sequence:
            self.build_graph()


    def build_graph(self):
        """-----------------------------------------------------------------------------------------
        Construct the De Bruijn graph with the specified kmer
        :return:
        -----------------------------------------------------------------------------------------"""
        k = self.k
        sequence = self.sequence
        nodes = self.nodes
        adj = self.adj
        for i in range(len(sequence) - k + 1):
            u = sequence[i: i + k - 1]
            v = sequence[i + 1: i + k]  # Actually next node is sequence[i+1 : i+k]
            adj[u].append(v)
            nodes.add(u)
            nodes.add(v)

        self.root = u
        return len(nodes)

    def arborescence(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        nodes = self.nodes
        last_exit = self.last_exit

        adj = self.adj

        # last_exit = {}
        in_tree = {self.root}

        for start_node in nodes:
            if start_node in in_tree: continue

            curr = start_node
            path = {}
            while curr not in in_tree:
                # Pick a random neighbor from the original transition list
                # Note: For exact k-mer counts, pick from available edges
                next_node = random.choice(adj[curr])
                path[curr] = next_node
                curr = next_node

            # Add the loop-erased path to the tree
            curr = start_node
            while curr not in in_tree:
                last_exit[curr] = path[curr]
                in_tree.add(curr)
                curr = path[curr]

        return

    def shuffle(self):
        """-----------------------------------------------------------------------------------------
        shuffle all edges EXCEPT the 'last-exit' edges
        :return:
        -----------------------------------------------------------------------------------------"""
        adj = self.adj
        last_exit = self.last_exit

        for u in adj:
            edges = adj[u]
            if u in last_exit:
                # Move the last_exit edge to the end of the list
                lastpos = edges.index(last_exit[u])
                edges[lastpos], edges[-1] = edges[-1],edges[lastpos]
                random.shuffle(edges[:-1])
                # edges.remove(last_exit[u])
                # random.shuffle(edges)
                # edges.append(last_exit[u])
            else:
                random.shuffle(edges)

            adj[u] = edges

        return

    def traverse(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        adj = self .adj
        curr = self.sequence[:k-1]
        result = [curr]

        while adj[curr]:
            next_node = adj[curr][0]
            adj[curr].remove(next_node)
            result.append(next_node[-1])  # Append the last character
            curr = next_node

        return "".join(result)
# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    k = 3
    alpha = 'ACGT'
    # sequence = ('AAAACCCAAAA')
    sequence = ('AAAACCCAAAA')
    g = DBG(sequence, 3)
    g.arborescence()
    g.shuffle()
    shuffled = g.traverse()
    print(f'{sequence} -> {shuffled}')

    exit(0)


    # for l in range(len):
    #     base = random.choice(alpha)
    #     sequence += base
    #     print(f'{base}:{sequence}')

    # overlap = defaultdict(Word)
    # words = defaultdict(Word)
    # prev = ''
    # for pos in range(0, len - k + 1):
    #     kmer = sequence[pos:pos + k]
    #     print(kmer)
    #
    #     words[kmer].label = kmer
    #     # overlap[kmer[1:-1]].out.append(words[kmer])
    #     if prev:
    #         words[prev].out.append(words[kmer])
    #         # overlap[prev[1:-1]].out.append(words[kmer])

    exit(0)

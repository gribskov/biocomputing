class Tree:
    """=============================================================================================
    A multifurcating tree class

    Michael Gribskov     15 April 2021
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        The tree is a list of linked nodes
        -----------------------------------------------------------------------------------------"""
        self.id = ''
        self.distance = None
        self.children = []

    def __iter__(self):
        return self.dfs()

    def dfs(self):
        """-----------------------------------------------------------------------------------------
        Generator function for depth-first traversal of tree
        :return: True
        -----------------------------------------------------------------------------------------"""
        yield self
        for child in self.children:
            yield from child.dfs()

        return True

    def balance(self, reverse=True):
        """-----------------------------------------------------------------------------------------
        Rearrange tree with the largest branches on the left (reverse=True) or on the right
        (reverse=False)
        :param reverse: boolean
        :return: True
        -----------------------------------------------------------------------------------------"""
        for node in self:
            if node.children:
                node.children.sort(key=lambda n: n.size(), reverse=reverse)

        return True

    def dump(self):
        """-----------------------------------------------------------------------------------------
        Dump the contents of the tree. Mostly for debugging

        :return: True
        -----------------------------------------------------------------------------------------"""

        print('\nnode:{}\n\tid:{}'.format(id(self), self.id))
        print('\tchildren:', end='')

        for child in self.children:
            print(' {}'.format(id(child)), end='')
        print()
        for child in self.children:
            child.dump()

        return True

    @staticmethod
    def find_comma(newick):
        """-----------------------------------------------------------------------------------------
        Find the central comma in the newick string
        this version only works for binary trees

        :param newick:
        :return: int, pos of comma in string, None if not found
        -----------------------------------------------------------------------------------------"""
        depth = 0
        pos = 0
        for c in newick:
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            elif c == ',':
                if depth == 0:
                    return pos

            pos += 1

        return None

    def from_newick(self, newick):
        """-----------------------------------------------------------------------------------------
        Construct the tree by recursively splitting the newick string at central commas until all
        commas have been removed
        :param newick: string
        :return: True
        -----------------------------------------------------------------------------------------"""
        branch_length, newick = self.get_branch_length(newick)
        if branch_length:
            self.id += str(branch_length)
        subtree = Tree.split_newick(newick)

        if subtree:
            for node in subtree:
                child = Tree()
                child.from_newick(node)
                self.children.append(child)

        else:
            # leaf node add string to id
            self.id = newick

        return True

    @staticmethod
    def get_branch_length(newick):
        """-----------------------------------------------------------------------------------------
        Check if a branch length is present, and if it is extract it, return the branch length
        and the shortened string.  if no branch length is present, return None and the original
        string

        :param newick: string
        :return: float, string
        -----------------------------------------------------------------------------------------"""
        last_colon = newick.rfind(':')
        last_paren = newick.rfind(')')
        if last_paren == -1 or last_paren > last_colon:
            return None, newick
        else:
            branch_length = float(newick[last_colon + 1:])
            newick = newick[:last_colon]
            return branch_length, newick

    def newick_string(self):
        """-----------------------------------------------------------------------------------------
        Recursively generate newick string for the tree
        :return:string
        -----------------------------------------------------------------------------------------"""
        newick = ''
        if self.children:
            newick += '('
            for child in self.children:
                newick += child.newick_string() + ','
            newick = newick.rstrip(',') + ')'
            if self.id:
                newick += ':' + self.id
        else:
            newick += self.id

        return newick

    def size(self):
        """-----------------------------------------------------------------------------------------
        Return number of nodes (including self)
        :return:
        -----------------------------------------------------------------------------------------"""
        n = 0
        for node in self:
            n += 1
        return n

    @staticmethod
    def split_newick(newick):
        """-----------------------------------------------------------------------------------------
        Split a newick string into its child nodes

        :param newick: string
        :return: list of newick substrings (child nodes of tree)
        -----------------------------------------------------------------------------------------"""
        # remove outer parentheses
        newick = newick.strip(' ;')
        if newick.startswith('('):
            newick = newick[1:-1]

        subtree = []
        comma = []

        depth = 0
        pos = 0
        for c in newick:
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            elif c == ',':
                if depth == 0:
                    comma.append(pos)
            pos += 1

        if comma:
            begin = 0
            for c in comma:
                subtree.append(newick[begin:c])
                begin = c + 1
            subtree.append(newick[begin:])

        return subtree


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # # test 1: does find_comma() work?
    # newick = [ '(a,b),(c, (d,e)', 'a,b', 'c, (d,e)']
    # for n in newick:
    #     comma = Tree.find_comma(n)
    #     print('{} => {}\t\t{}'.format(n, n[:comma], n[comma+1:]))

    # test 2: build a tree
    # newick = ['((a,b),(c, (d,e)));']
    # for n in newick:
    #     tree = Tree()
    #     tree.from_newick(n)
    #     tree.dump()

    newick = ['((a,b),(c, (d,e)));',
              '((raccoon,bear),((sea_lion, seal),((monkey,cat), weasel)),dog);',
              '((raccoon:19.2,bear:6.8):0.9,((sea_lion:13.0, seal:12.0):7.5,((monkey:100.9,cat:47.1):20.6, weasel:18.9):2.09):3.9,dog:25.5);',
              '(dog:25.5,(raccoon:19.2,bear:6.8):0.9,((sea_lion:13.0, seal:12.0):7.5,((monkey:100.9,cat:47.1):20.6,weasel:18.9):2.09):3.9);']
    for n in newick:
        tree = Tree()
        tree.from_newick(n)
        # tree.dump()

        print('tree size: {}'.format(tree.size()))
        for node in tree:
            print('node:{}\tid:{}:'.format(id(node), node.id))

        print(tree.newick_string())

        tree.balance()
        print(tree.newick_string())
        tree.balance(reverse=False)
        print(tree.newick_string())

    exit(0)

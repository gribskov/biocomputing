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
        # remove outer parentheses
        newick = newick.strip(' ;')
        if newick.startswith('('):
            newick = newick[1:-1]

        # find central comma
        comma = Tree.find_comma(newick)

        if comma:
            # build child nodes
            left = newick[:comma]
            right = newick[comma + 1:]

            child = Tree()
            child.from_newick(left)
            self.children.append(child)

            child = Tree()
            child.from_newick(right)
            self.children.append(child)
        else:
            # leaf node add string to id
            self.id = newick

        return True


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
    newick = ['((a,b),(c, (d,e)));']
    for n in newick:
        tree = Tree()
        tree.from_newick(n)
        tree.dump()

    exit(0)

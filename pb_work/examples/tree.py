"""=================================================================================================
Tree class: 2019 practical biocomputing example

Michael Gribskov     16 February 2019
================================================================================================="""


class Tree:
    """=============================================================================================
    Represents binary trees with branch lengths.  The object is actually just a node in the tree.

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        A tree node has the following attributes
            name (string)
            branch_length (float)
            left (Tree) - the left child
            right (Tree) - the right child
        -----------------------------------------------------------------------------------------"""
        self.name = ''
        self.branch_length = 0.0
        self.left = None
        self.right = None

    def dump(self):
        print('\n{}'.format(self))
        print('\tname:{}'.format(self.name))
        if self.branch_length:
            print('\tbranch length:{}'.format(self.branch_length))
        print('\tleft:{}'.format(self.left))
        print('\tright:{}'.format(self.right))

        # recursively dump child nodes
        if self.left:
            self.left.dump()

        if self.right:
            self.right.dump()

    def tree_length(self):
        """-----------------------------------------------------------------------------------------
        Recursively calculate the length of the tree rooted at the  current node
        :return: float (total length of tree)
        -----------------------------------------------------------------------------------------"""
        l = self.branch_length
        if self.left:
            l += self.left.tree_length()

        if self.right:
            l += self.right.tree_length()

        return l

    def nodes(self):
        """-----------------------------------------------------------------------------------------
        Return a list of nodes in depth-first order, rooted at this node

        :return: list of Tree
        -----------------------------------------------------------------------------------------"""
        nodelist = [self]
        if self.left:
            nodelist += self.left.nodes()

        if self.right:
            nodelist += self.right.nodes()

        return nodelist


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    root = Tree()

    print('make a single node named root')
    root.name = 'root'
    root.dump()

    print('\nadd children ab and cde')
    ab = Tree()
    ab.name = 'ab'
    ab.branch_length = 1.1

    cde = Tree()
    cde.name = 'cde'
    cde.branch_length = 0.5

    root.left = cde
    root.right = ab

    print('\nadd children of cde')
    c = Tree()
    c.name = 'c'
    c.branch_length = 0.7

    de = Tree()
    de.name = 'de'
    de.branch_length = 0.5

    cde.right = c
    cde.left = de

    root.dump()
    print('\ntree length:', root.tree_length())

    print('\nNode list')
    print('root:')
    l = root.nodes()
    for node in root.nodes():
        print('\t{}'.format(node.name))

    print('cde:')
    for node in cde.nodes():
        print('\t{}'.format(node.name))

    exit(0)

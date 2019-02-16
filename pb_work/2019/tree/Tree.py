"""=================================================================================================
Tree class

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

        return True


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    root = Tree()

    print('make a single node named root')
    root.name = 'root'
    root.dump()

    print('add children ab and cde')
    ab = Tree()
    ab.name = 'ab'
    ab.branch_length = 1.1

    cde = Tree()
    cde.name = 'cde'
    cde.branch_length = 0.5

    root.left = cde
    root.right = ab

    root.dump()
    root.left.dump()
    root.right.dump()


    exit(0)

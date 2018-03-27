class Tree():
    """=============================================================================================
    Class for multifurcating trees.

    27 March 2018     Michael Gribskov
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Tree constructor.
        -----------------------------------------------------------------------------------------"""
        self.name = ''
        self.children = []

    def __str__(self):
        """-----------------------------------------------------------------------------------------
        return a formatted string describing a tree node
        :return: string
        -----------------------------------------------------------------------------------------"""
        treestr = 'node:{}\n    children:'.format(self.name)
        if self.children:
            treestr += ','.join(str(self.children))

        else:
            treestr += 'None'

        return treestr


# ==================================================================================================
# testing
# ==================================================================================================
if __name__ == '__main__':
    # create a tree
    tree = Tree()

    # test entry of information into tree node
    tree.name = 'a'
    print(str(tree))

    # add child nodes to a
    b = Tree()
    b.name = 'b'
    c = Tree()
    c.name = 'c'

    tree.children.append(b)
    tree.children.append(c)
    print(str(tree))

    exit(0)

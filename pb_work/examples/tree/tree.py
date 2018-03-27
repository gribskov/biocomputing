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
        treestr = 'node:{}\n    children: '.format(self.name)
        if self.children:
            for child in self.children:
                treestr += str(child.name) + ' '

        else:
            treestr += 'None'

        return treestr

    def dfs(self):
        """-----------------------------------------------------------------------------------------
        generator to trace the tree in depth first order
        :yield: next tree node
        -----------------------------------------------------------------------------------------"""
        yield self
        for child in self.children:
            for node in child.dfs():
                yield node

        StopIteration

# ==================================================================================================
# testing
# ==================================================================================================
if __name__ == '__main__':
    # create a tree
    tree = Tree()

    # test entry of information into tree node
    print('\n-------------\ntree a with children b, c\n-------------')

    tree.name = 'a'
    print(str(tree))

    # add child nodes to a
    b = Tree()
    b.name = 'b'
    c = Tree()
    c.name = 'c'

    tree.children.append(b)
    tree.children.append(c)
    print('\n{}'.format(str(tree)))

    # test dfs iterator
    print('\n-------------\ndfs iteration\n-------------')
    for node in tree.dfs():
        print('\n{}'.format(str(node)))

    exit(0)

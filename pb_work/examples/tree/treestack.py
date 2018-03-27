"""=================================================================================================
Stack based approach to tree iteration

27 march 2018     Michael Gribskov
================================================================================================="""
from tree import Tree


def dfsStack(tree, fxn):
    """---------------------------------------------------------------------------------------------
    Depth first traversal of tree.
    :param tree:
    :param txn: a function to execute at each node of the tree
    :return: integer, number of nodes visited
    ---------------------------------------------------------------------------------------------"""
    stack = []
    stack.append(tree)

    n = 0
    while stack:
        node = stack.pop()
        fxn(node)
        n += 1

        for child in reversed(node.children):
            stack.append(child)

    return n


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    # create a tree
    tree = Tree('a')

    # add child nodes to a
    b = tree.childAdd('b')
    c = tree.childAdd('c')
    print('\n{}'.format(str(tree)))

    # add a child to b
    b.childAdd('d')

    # test dfs iterator
    print('\n-------------\ndfs iteration (main)\n-------------')
    print('expected order:a, b, d, c')

    stack = []
    stack.append(tree)

    while stack:
        node = stack.pop()
        print('node:{}'.format(node.name))

        for child in reversed(node.children):
            stack.append(child)

    print('\n-------------\ndfs iteration (function)\n-------------')
    print('expected order:a, b, d, c')


    def printnode(node):
        print('node {}'.format(node.name))


    dfsStack(tree, printnode)

    print('\nusing lambda expression')
    dfsStack(tree, lambda n: print('node {}'.format(n.name)))

    exit(0)

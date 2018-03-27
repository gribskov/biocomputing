class Tree():
    """=============================================================================================
    Class for multifurcating trees.

    27 March 2018     Michael Gribskov
    ============================================================================================="""
    import sys

    def __init__(self, name=''):
        """-----------------------------------------------------------------------------------------
        Tree constructor.
        -----------------------------------------------------------------------------------------"""
        self.name = name
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

    def childAdd(self, name=''):
        """-----------------------------------------------------------------------------------------
        Create a new child node AND add to children of current node.  Return the new node
        :return:
        -----------------------------------------------------------------------------------------"""
        new = Tree(name)
        self.children.append(new)

        return new

    def newickLoad(self, newick):
        """-----------------------------------------------------------------------------------------
        Builds a tree from a newick string.  Here are some examples from Joe Felsenstein's website
        http://evolution.genetics.washington.edu/phylip/newicktree.html. These are the test trees
        used below

        unrooted tree (three central branches)
        ((raccoon:19.2,bear:6.8):0.9,((sea_lion:13.0, seal:12.0):7.5,((monkey:100.9,cat:47.1):20.6,
        weasel:18.9):2.09):3.9,dog:25.5);

        rooted tree
        (Bovine:0.7,(Gibbon:0.4,(Orang:0.3,(Gorilla:0.2,(Chimp:0.2, Human:0.1):0.1):0.1):0.2):0.5,
        Mouse:1.2):0.1;

        another unrooted tree
        (Bovine:0.69,(Hylobates:0.36,(Pongo:0.34,(G._Gorilla:0.17, (P._paniscus:0.19,
        H._sapiens:0.12):0.08):0.06):0.15):0.55, Rodent:1.21);

        USAGE
           root.newickLoad( newick_string );

        :param newick: a newick formatted string defining a tree
        :return: True, if successful
        -----------------------------------------------------------------------------------------"""
        newick = newick.strip()  # removes white space
        newick = newick.rstrip(';')  # removes trailing semicolon if present
        if newick[0] != '(':
            sys.stderr.write('Tree.newickLoad - Newick string must begin with (\n')
            return False

        # push the root node on the stack
        stack = [self]
        node = self
        word = ''

        for letter in newick:
            if letter.isspace():
                continue
            elif letter == '(':
                newnode = node.childAdd()
                stack.append(node)
                node = newnode
                word = ''

            elif letter == ',':
                node.name += word
                node = stack[-1]
                newnode = node.childAdd()
                node = newnode
                word = ''

            elif letter == ')':
                node.name += word
                node = stack.pop()
                word = ''

            else:
                word += letter

        # when you finish, if there is anything in word it belongs to the root
        node.name += word

        return True

    # end of newickLoad

    def newick(self):
        """-----------------------------------------------------------------------------------------
        generate newick string
        newick format show the tree as a set of nested parentheses in which the children of a node
        are shown as a comma delimited list inside a pair of parenthess.
        a --|
            |--|
        b --|  |
               |--
        c -----|
        ((a,b),c)

        :return: newick formatted tree
        -----------------------------------------------------------------------------------------"""
        newick = ''
        punct = '('
        if self.children:
            for child in self.children:
                newick += punct + child.newick()
                punct = ','
            newick += '){}'.format(self.name)
        else:
            newick = self.name

        return newick

# ==================================================================================================
# testing
# ==================================================================================================
if __name__ == '__main__':

    # test entry of information into tree node
    print('\n-------------------------\ntree a with children b, c\n-------------------------')

    # create a tree
    tree = Tree('a')

    # add child nodes to a
    b = tree.childAdd('b')
    c = tree.childAdd('c')
    print('\n{}'.format(str(tree)))

    # add a child to b
    b.childAdd('d')

    # test dfs iterator
    print('\n-------------\ndfs iteration\n-------------')
    print('expected order:a, b, d, c')

    for node in tree.dfs():
        print('\n{}'.format(str(node)))

    # read newick formatted tree
    tree = ' ((a,b),c);'
    title = 'newick formatted tree' + tree
    print('\n{}\n{}\n{}'.format('-' * len(title), title, '-' * len(title)))
    n1 = Tree()
    n1.newickLoad(tree)
    for node in n1.dfs():
        print('\n{}'.format(str(node)))

    tree = ' (Bov:0.69,(Hyl:0.36,(Pon:0.34,(Gor:0.17,(Pan:0.19,Hom:0.12):0.08):0.06):0.15):0.55,Rod:1.21);'
    title = 'newick formatted tree' + tree
    print('\n{}\n{}\n{}'.format('-' * len(title), title, '-' * len(title)))
    n2 = Tree()
    n2.newickLoad(tree)
    for node in n2.dfs():
        print('\n{}'.format(str(node)))

    title = 'write newick formatted tree' + tree
    print('\n{}\n{}\n{}'.format('-' * len(title), title, '-' * len(title)))

    print(n2.newick())

    exit(0)

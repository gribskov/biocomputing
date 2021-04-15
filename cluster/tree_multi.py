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

    def from_newick(self, newick):
        """-----------------------------------------------------------------------------------------
        Construct the tree by recursively splitting the newick string at central commas until all
        commas have been removed
        :param newick: string
        :return: True
        -----------------------------------------------------------------------------------------"""
        # remove outer parentheses

        # find central comma

        # build child nodes

        return True


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    exit(0)

"""=================================================================================================
Two classes in this file.
    GeneOntologyItem: parses an OBO formatted entry for a single GO term and stores is
    GeneOntology: a container of GeneOntologyItems that can iterate and load an entire GO file

Michael Gribskov     12 January 2019
================================================================================================="""
import sys


class GeneOntologyItem():
    """=============================================================================================
    the object is an array of information tag-value pairs, stored in GeneOntologyItem.info
    ============================================================================================="""

    def __init__(self, obo_block):
        """-----------------------------------------------------------------------------------------
        GeneOntologyItem constructor
        -----------------------------------------------------------------------------------------"""
        self.type = ''
        self.info = {}
        self.count = {}
        if obo_block:
            self.parse(obo_block)

    def parse(self, text):
        """-----------------------------------------------------------------------------------------
        parse a block of OBO formatted text, stored as a list of strings.
        each line is a tag value pair separated by the first colon

        the keys are used as dictionary keys and the values are stored as lists of strings

        :param text: list of lines read from file
        :return: number of liens processed
        -----------------------------------------------------------------------------------------"""
        nread = 0

        # the first line determines the record type
        self.type = 'Header'
        if '[' in text[0]:
            self.type = text[0].lstrip('[').rstrip(']\n')

        for pair in text[1:]:
            nread += 1
            key, value = pair.strip(' \n').split(':', maxsplit=1)
            key = key.strip(' ')
            value = value.strip(' ')

            if key in self.info:
                self.info[key].append(value)
                self.count[key] += 1
            else:
                self.info[key] = [value]
                self.count[key] = 1

        self.n = nread

        return nread

    def id(self):
        """-----------------------------------------------------------------------------------------
        REturn the GO ID if present.  This should be in a line such as
        id: GO:0000006
        in the OBO file, and should be stored in the info as a single key value pair

        :return: string, empty string indicates not present
        -----------------------------------------------------------------------------------------"""
        id = ''
        if 'id' in self.info:
            if len(self.info['id']) > 1:
                sys.stderr.write('GeneOntologyItem::id - ID multiply defined {}'.
                                 format(self.info['id']))
            id = self.info['id'][0]

        return id

    def feature(self, key):
        """-----------------------------------------------------------------------------------------
        Return the array of values for any defined feature.  If feature is unknown, and empty
        list is returned.  An empty list tests as False

        :param key: string
        :return: list of strings or empty list
        -----------------------------------------------------------------------------------------"""
        result = ()
        if key in self.info:
            result = self.info[key]

        return result


class GeneOntology(object):
    """=============================================================================================
    Gene ontology container
    ============================================================================================="""

    def __init__(self, file=''):
        """-----------------------------------------------------------------------------------------
        GeneOntology is mostly a bag of GeneOntologyItem. The list of GO terms are are stored in
        term.
        index is a dictionary indexed by the GO ID, e.g. GO:0051082, that allows you to navigate
        between terms following relationships such as is_a
        -----------------------------------------------------------------------------------------"""
        self.filename = file
        self.fh = None
        self.block = []
        self.term = []
        self.index = {}
        self.parent = {}
        self.child = {}

        if self.filename:
            if not self.openfile():
                sys.stderr.write('\nGeneOntology:init - unable to open file ({})\n'.
                                 format(self.filename))

    def openfile(self):
        """-----------------------------------------------------------------------------------------
        Safe file open
        if successful returns True and sets the filehandle
        if unsuccessful returns false, but does not die

        :return: logical
        -----------------------------------------------------------------------------------------"""
        try:
            self.fh = open(self.filename, 'r')
        except (OSError, IOError):
            return False

        return True

    def read_obo_block(self):
        """-----------------------------------------------------------------------------------------
        OBO format is broken up into blocks separated by blank lines.  This function reads the
        block beginning at the current file position. Each block corresponds to a GeneOntologyItem.

        :return: logical, True if block is read
        "----------------------------------------------------------------------------------------"""
        status = False
        self.block = []

        line = self.fh.readline()
        if not line:
            # end of file
            return status

        while line.rstrip(' \n'):
            status = True
            self.block.append(line)
            line = self.fh.readline()

        return True

    def next(self):
        """-----------------------------------------------------------------------------------------
        parse the current block and read the next one from the open filehandle
        Returns False when no block is read

        :return: logical
        -----------------------------------------------------------------------------------------"""
        block_available = self.read_obo_block()
        if block_available:
            item = GeneOntologyItem(self.block)
            self.term.append(item)
            id = item.id()
            idx = len(self.term) - 1
            if id:
                # index GO terms
                self.index[id] = idx
            alt_id = item.feature('alt_id')
            if alt_id:
                for id in alt_id:
                    self.index[id] = idx

        return block_available

    def load(self):
        """-----------------------------------------------------------------------------------------
        Read the entire obo file and add to the object.  does not clear internal storage so multiple
        files can be read into the same object

        :return: int, number of records loaded (all types)
        -----------------------------------------------------------------------------------------"""
        n_load = 0
        while self.next():
            n_load += 1
            if not n_load % 1000:
                print('.', end='')

        return n_load

    def parent_child(self):
        """-----------------------------------------------------------------------------------------
        Identify parent child relationships, and build two indexes base on the serial index

        :return: it, number of parent-child relationships
        -----------------------------------------------------------------------------------------"""
        idx = 0
        for term in self.term:
            if not term.type == 'Term':
                idx += 1
                continue

            print('index: {}\tid: {}\tname {}'.format(
                idx, term.info['id'][0], term.info['name'][0]))
            idx += 1
            if 'is_a' not in term.info:
                continue

            for is_a in term.info['is_a']:
                go_id, comment = is_a.split('!', maxsplit=1)
                go_id = go_id.strip(' ')
                comment = comment.strip(' ')
                print('\tgo_id: {}\tcomment: {}'.format(go_id, comment))
                try:
                    self.parent[term.info['id'][0]].append(go_id)
                except KeyError:
                    self.parent[term.info['id'][0]] = [go_id]

                try:
                    self.child[go_id].append(term.info['id'][0])
                except KeyError:
                    self.child[go_id] = [term.info['id'][0]]

        return True

    def trace(self,leaf):
        """-----------------------------------------------------------------------------------------
        trace the path from a leaf to a node with no parent (i.e., the root)

        :param leaf: string, a GO ID
        :return: list of string, GO IDs
        -----------------------------------------------------------------------------------------"""
        golist = []
        while leaf in self.parent:
            golist.append(leaf)
            idx = self.index[leaf]
            term = self.term[idx]
            print('{} {} ({})'.format(term.info['id'][0], term.info['name'][0], idx))

            if leaf in self.parent:
                leaf = self.parent[leaf][0]

        golist.append(leaf)
        idx = self.index[leaf]
        term = self.term[idx]
        print('{} {} ({})'.format(term.info['id'][0], term.info['name'][0], idx))

        return True

# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gofile = 'go-basic.obo'
    go = GeneOntology(file=gofile)
    # go.next()
    print('file: {}'.format(go.filename))
    nloaded = go.load()
    print('{} terms loaded from {}'.format(nloaded, gofile))

    all_features = {}
    for goitem in go.term:
        for key in goitem.info:
            if key in all_features:
                all_features[key] += 1
            else:
                all_features[key] = 1

    print('\nFeatures')
    for key in all_features:
        print('{}: {}'.format(key, all_features[key]))

    go.parent_child()

    n = 0
    for c in go.child:
        if c not in go.parent:
            print('index: {}'.format(go.index[c]))
            term = go.term[go.index[c]]
            if 'obsolete' not in term.info:
                n += 1
                print('{} child but no parent: {}'.format(n, c))
                print(term.info)

    gotest = ['GO:2001314', 'GO:0014041', 'GO:1902679']
    for g in gotest:
        print('\nstart:{}'.format(g))
        go.trace(g)

    exit(0)

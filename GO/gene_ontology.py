"""=================================================================================================


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
            key, value = pair.rstrip('\n').split(':', maxsplit=1)

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
        GeneOntology is mostly a bag of GeneOntologyItem. The GO terms are are stored in term.
        index is a dictionary in dexed by the GO ID, e.g. GO:0051082, that allows you to navigate
        between terms following relationships such as is_a
        -----------------------------------------------------------------------------------------"""
        self.filename = file
        self.fh = None
        self.block = []
        self.term = []
        self.index = {}

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
            if id:
                # index GO terms
                self.index[id] = len(self.term)

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

    exit(0)

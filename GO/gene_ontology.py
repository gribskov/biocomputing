"""=================================================================================================


Michael Gribskov     12 January 2019
================================================================================================="""
import sys


class GeneOntologyItem():
    """=============================================================================================
    the object is an array of information tag-value pairs, stored in GeneOntologyItem.info
    ============================================================================================="""

    def __init__(self, GO):
        """-----------------------------------------------------------------------------------------
        GeneOntologyItem constructor
        -----------------------------------------------------------------------------------------"""
        self.type = ''
        self.info = []
        self.count = {}
        if GO.block:
            self.parse(GO.block)


    def parse(self, text):
        """-----------------------------------------------------------------------------------------
        parse a block of OBO formatted text, stored as a list of strings.
        each line is a tag value pair separated by the first colon

        returns the key-value pairs as a list of ditionaries in self.info
        returns the count of each key in self.count

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
            self.info.append({'key': key, 'value': value})
            try:
                self.count[key] += 1
            except KeyError:
                self.count[key] = 1

        self.n = nread

        return nread


class GeneOntology(object):
    """=============================================================================================
    Gene ontology container
    ============================================================================================="""

    def __init__(self, file=''):
        """-----------------------------------------------------------------------------------------
        GeneOntology is mostly a bag of GeneOntologyItem
        -----------------------------------------------------------------------------------------"""
        self.filename = file
        self.fh = None
        self.block = []
        self.term = []

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
            self.obo_block()
        except (OSError, IOError):
            return False

        return True

    def obo_block(self):
        """-----------------------------------------------------------------------------------------
        OBO format is broken up into blocks separated by blank lines.  This function reads the
        block beginning at the current file position

        :return: logical, True if block is read
        "----------------------------------------------------------------------------------------"""
        status = False
        self.block = []

        line = self.fh.readline()
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
        item = GeneOntologyItem(self)
        self.term.append(item)
        while not item.type == 'Term':
            if self.obo_block():
                item = GeneOntologyItem(self)
                self.term.append(item)
            else:
                break

        return self.obo_block()


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gofile = 'go-basic.obo'
    go = GeneOntology(file=gofile)
    go.next()
    print('file: {}'.format(go.filename))

    exit(0)

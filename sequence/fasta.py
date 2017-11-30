'''
Fasta sequence class.  Supports iteration over a multi-fasta file
    id
    documentation
    sequence
'''


# import re

class Fasta:
    def __init__(self):
        self.filename = ''
        self.id = ''
        self.doc = ''
        self.seq = ''
        self.buffer = ''

    def open(self, filename):
        '''
        open a file for reading
        '''
        #print('opening:', filename)
        try:
            self.fh = open(filename, 'r')
        except:
            pass

        self.filename = filename

    def next(self):
        '''
        return the next entry from an open file into the object
        '''
        self.read();
        if not self.id:
            return 0
        else:
            return 1

    def read(self):
        '''
        read one sequence from the file, leave the following line in buffer
        usage:
        fasta.read()
        '''

        self.id = ''
        self.doc = ''
        self.seq = ''

        self.getID()

        for line in self.fh:
            if line.isspace(): continue
            line = line.rstrip('\n')

            if line[0] == '>':
                self.buffer = line
                return
            else:
                self.seq += line

    def getID(self):
        '''
        intended to be used internally for sequence reading
        check the buffer, and if not empty read the ID and documentation
        sore in id and doc attributes
        id will be stripped of >
        documentation will be and empty string if there is nothing following the ID
        '''
        # if buffer is empty read a line
        if self.buffer:
            line = self.buffer
            self.buffer = ''
        else:
            line = self.fh.readline()
            line = line.rstrip('\n')

        # get the ID and documentation from the doc line
        try:
            id, doc = line.split(" ", 1)
        except ValueError:
            # documentation is missing
            id = line
            doc = ''

        self.id = id.lstrip('>')
        self.doc = doc

    def length(self):
        '''
        return the length of the  current sequence
        return 0 if there is none
        usage
            seqlen = fasta.length()
        '''
        return len(self.seq)

    def format(self):
        '''
        return a formatted string with the current sequence
        usage
            seq = fasta.format()
        '''
        # hardwired length for sequence lines
        # TODO change to argument
        LINELEN = 50

        string = '>{0} {1}\n'.format(self.id, self.doc)
        pos = 0
        while pos < len(self.seq):
            string += '{0}\n'.format(self.seq[pos:pos + LINELEN])
            pos += LINELEN

        return string

    def trimDocByRegex(self, target):
        '''
        Shorten documentation by substituting the target regex with nothing
        target must be a compiled regex
        The new documentation string is returned
        usage
            trim = re.compile( 'len=\d+ ' )
            doc = fasta.trimDocAfterMatch( trim )
        '''
        self.doc = target.sub('', self.doc)

        return self.doc

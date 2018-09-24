'''-------------------------------------------------------------------------------------------------
Fasta sequence class.  Supports iteration over a multi-fasta file
    filename
    id
    documentation
    sequence

    Synopsis
    from sequence.fasta import Fasta

    fasta = Fasta()
    fasta.open('filename')
    while fasta.next():
        print(fasta.format(linelen=60))

-------------------------------------------------------------------------------------------------'''


class Fasta:

    def __init__(self, file=""):
        """-----------------------------------------------------------------------------------------
        Fasta class constructor. Attributes
            filename
            id
            doc
            seq
            buffer  (read ahead buffer, only internal)
        -----------------------------------------------------------------------------------------"""
        self.filename = file
        self.id = ''
        self.doc = ''
        self.seq = ''
        self.buffer = ''

        if self.filename:
            self.open(self.filename)

    def open(self, filename):
        '''-----------------------------------------------------------------------------------------
        open a file for reading
        -----------------------------------------------------------------------------------------'''
        # print('opening:', filename)
        try:
            self.fh = open(filename, 'r')
        except:
            print('Fasta::open - file open error')

        self.filename = filename

    def next(self):
        '''-----------------------------------------------------------------------------------------
        return the next entry from an open file into the object
        usage
            while fasta.next():
               ...
        -----------------------------------------------------------------------------------------'''
        return self.read()

    def read(self):
        '''-----------------------------------------------------------------------------------------
        read one sequence from the file, leave the following line in buffer
        usage:
        fasta.read()
        -----------------------------------------------------------------------------------------'''

        self.id = ''
        self.doc = ''
        self.seq = ''

        # get the ID and doc
        self.getID()

        for line in self.fh:
            if line.isspace(): continue
            try:
                line = line.rstrip('\n')
            except TypeError:
                # in case of a byte string
                line = line.decode()
                line = line.rstrip('\n')

            if line[0] == '>':
                self.buffer = line
                break

            else:
                self.seq += line

        if len(self.id) > 0 or len(self.doc) > 0 or len(self.seq) > 0:
            return True

        # fall through to false if nothing can be read
        return False

    def getID(self):
        '''-----------------------------------------------------------------------------------------
        intended to be used internally for sequence reading
        check the buffer, and if not empty read the ID and documentation
        sore in id and doc attributes
        id will be stripped of >
        documentation will be and empty string if there is nothing following the ID
        -----------------------------------------------------------------------------------------'''
        # if buffer is empty read a line
        if self.buffer:
            line = self.buffer
            self.buffer = ''
        else:
            line = self.fh.readline()
            try:
                line = line.rstrip('\n')
            except TypeError:
                # in case of a byte string
                # line = line.rstrip('\n'.encode())
                line = line.decode()
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
        '''-----------------------------------------------------------------------------------------
        return the length of the  current sequence
        return 0 if there is none
        usage
            seqlen = fasta.length()
        -----------------------------------------------------------------------------------------'''
        return len(self.seq)

    def format(self, linelen=50):
        '''-----------------------------------------------------------------------------------------
        return a formatted string with the current sequence
        usage
            seq = fasta.format()
        -----------------------------------------------------------------------------------------'''
        string = '>{0} {1}\n'.format(self.id, self.doc)
        pos = 0
        while pos < len(self.seq):
            string += '{0}\n'.format(self.seq[pos:pos + linelen])
            pos += linelen

        return string

    def trimDocByRegex(self, target):
        '''-----------------------------------------------------------------------------------------
        Shorten documentation by substituting the target regex with nothing
        target must be a compiled regex
        The new documentation string is returned
        usage
            trim = re.compile( 'len=\d+ ' )
            doc = fasta.trimDocAfterMatch( trim )
        -----------------------------------------------------------------------------------------'''
        self.doc = target.sub('', self.doc)

        return self.doc

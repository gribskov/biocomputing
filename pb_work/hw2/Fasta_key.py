"""=================================================================================================
Fasta sequence class.  Supports iteration over a multi-fasta file
    filename
    id
    documentation
    sequence

    Synopsis
    from fasta import Fasta

    fasta = Fasta()
    fasta.open('filename')
    while True:
        fasta.next():
        print(fasta.format(linelen=60))

================================================================================================="""


class Fasta:

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Fasta class constructor. Attributes
            filename
            id
            doc
            seq
            buffer  (read ahead buffer, only internal)
        -----------------------------------------------------------------------------------------"""
        self.filename = ''
        self.fh = None
        self.id = ''
        self.doc = ''
        self.seq = ''
        self.buffer = ''

        return

    def open(self, filename):
        """-----------------------------------------------------------------------------------------
        open a file for reading

        :return: True if successful, False otherwise
        -----------------------------------------------------------------------------------------"""
        status = True
        try:
            self.fh = open(filename, 'r')
        except FileNotFoundError:
            print('Fasta::open - file open error')
            status = False

        self.filename = filename

        return status

    def next(self):
        """-----------------------------------------------------------------------------------------
        read one sequence from the file, leave the following line in buffer
        usage:
            fasta.read()

        :return: True (successful), False (usually EOF)
        -----------------------------------------------------------------------------------------"""

        self.id = ''
        self.doc = ''
        self.seq = ''

        # get the ID and doc
        self.getID()

        for line in self.fh:
            if line.isspace():
                continue

            line = line.rstrip('\n')
            if line[0] == '>':
                self.buffer = line
                break

            else:
                line = line.replace('N','')
                line = line.replace('*','')
                self.seq += line

        if len(self.id) > 0 or len(self.doc) > 0 or len(self.seq) > 0:
            return True

        # fall through to false if nothing can be read
        return False

    def getID(self):
        """-----------------------------------------------------------------------------------------
        intended to be used internally for sequence reading
        check the buffer, and if not empty read the ID and documentation
        sore in id and doc attributes
        id will be stripped of >
        documentation will be and empty string if there is nothing following the ID

        :return: None
        -----------------------------------------------------------------------------------------"""
        # if buffer is empty read a line
        line = ''
        if self.buffer:
            line = self.buffer
            self.buffer = ''
        else:
            for line in self.fh:
                if line.isspace():
                    continue
                break
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

        return None

    def format(self, linelen=50):
        """-----------------------------------------------------------------------------------------
        return a formatted string with the current sequence
        usage
            seq = fasta.format()

        :param linelen: number of sequence characters per line
        -----------------------------------------------------------------------------------------"""
        string = '>{0} {1}\n'.format(self.id, self.doc)
        pos = 0
        while pos < len(self.seq):
            string += '{0}\n'.format(self.seq[pos:pos + linelen])
            pos += linelen

        return string


#  =================================================================================================
# main program - do not change the code below
# ==================================================================================================
if __name__ == '__main__':

    fasta = Fasta()

    fasta.open('short_broke.fa')
    while True:
        if fasta.next():
            print('sequence {}:{}\nlength={}'.format(fasta.id, fasta.doc, len(fasta.seq)))
            print('{}'.format(fasta.format(linelen=50)))
        else:
            break

    exit(0)

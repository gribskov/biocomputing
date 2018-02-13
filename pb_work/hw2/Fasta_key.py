"""=================================================================================================
Fasta sequence class.  Supports iteration over a multi-fasta file.  Object have the following
externally accessible attributes:

    filename    input filename
    fh          input filehandle
    id          ID string (without>)
    doc         documentation string
    seq         sequence (one letter code, sequence characters only)

    Synopsis
    from fasta import Fasta

    fasta = Fasta()
    fasta.open('filename')
    while True:
        fasta.next():
        print(fasta.format(linelen=60))

TODO: Don't forget to put your name and the date here
================================================================================================="""


class Fasta:

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Fasta class constructor. Attributes
            filename
            fh
            id
            doc
            seq
            _buffer  (read ahead buffer, internal use only)
        -----------------------------------------------------------------------------------------"""
        self.filename = ''
        self.fh = None
        self.id = ''
        self.doc = ''
        self.seq = ''
        self._buffer = ''

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
            print('Fasta::open - {} not found'.format(filename))
            status = False

        self.filename = filename

        return status

    def next(self):
        """-----------------------------------------------------------------------------------------
        read one sequence from the file, leave the following line in _buffer
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
                self._buffer = line
                break

            else:
                line = line.replace('N', '')
                line = line.replace('*', '')
                self.seq += line

        if len(self.id) > 0 or len(self.doc) > 0 or len(self.seq) > 0:
            return True

        # fall through to false if nothing can be read
        return False

    def getID(self):
        """-----------------------------------------------------------------------------------------
        intended to be used internally for sequence reading
        check the _buffer, and if not empty read the ID and documentation
        sore in id and doc attributes
        id will be stripped of >
        documentation will be and empty string if there is nothing following the ID

        :return: None
        -----------------------------------------------------------------------------------------"""
        # if _buffer is empty read a line
        line = ''
        if self._buffer:
            line = self._buffer
            self._buffer = ''
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
        :return: string containing formatted sequence
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

    outfile = 'hw2_key.out'
    try:
        out = open(outfile, 'w')
    except:
        print('unable to open {}'.format(outfile))
        exit(1)

    fasta = Fasta()
    fasta.open('short_broke.fa')

    n = 0
    while True:
        n += 1
        if fasta.next():
            print('sequence {}:{} length={}'.format(n, fasta.id, len(fasta.seq)))
            out.write('{}'.format(fasta.format(linelen=50)))
        else:
            break

        if n >= 10:
            break

    out.close()
    exit(0)

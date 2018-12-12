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
    codon2aa = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
                "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
                "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
                "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",

                "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
                "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
                "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
                "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",

                "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
                "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
                "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
                "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",

                "TAA": "_", "TAC": "Y", "TAG": "_", "TAT": "T",
                "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
                "TGA": "_", "TGC": "C", "TGG": "W", "TGT": "C",
                "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}

    complement = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a'}

    def __init__(self, file=""):
        """-----------------------------------------------------------------------------------------
        Fasta class constructor. Attributes
            filename
            id
            doc
            seq
            buffer  (read ahead buffer, only internal)
        -----------------------------------------------------------------------------------------"""
        self.filename = ''
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
        string = '>{0} {1}'.format(self.id, self.doc)
        pos = 0
        while pos < len(self.seq):
            string += '\n{0}'.format(self.seq[pos:pos + linelen])
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

    def reverseComplement(fasta):
        """-----------------------------------------------------------------------------------------
        Return the sequence converted to reverse complement
        :return: string
        -----------------------------------------------------------------------------------------"""
        seq = fasta.seq
        seq = seq.translate(Fasta.complement)

        return seq[::-1]

    def translate(fasta, frame=0, direction='f'):
        """-----------------------------------------------------------------------------------------
        translate in a nucleic acid sequence in the desired direction (f,r) and frame (0..2)
        incomplete codons at end are not translated
        stop codons are shown as '*'
        codons with ambiguity characters are translated as 'X';
        currently uses standard amino acid code
        TODO use supplied genetic code
        :param frame: integer, offset from beginning of sequence
        :param direction: string, forward (f) or reverse(r)
        :return: Fasta object (new)
        -----------------------------------------------------------------------------------------"""
        if not fasta.isACGT():
            sys.stderr.write('Fasta::translate - sequence must be ACGT')

        rf = '{}{}'.format(direction, frame)
        trans = Fasta()
        trans.id = fasta.id + '_{}'.format(rf)
        trans.doc = fasta.id + ' reading_frame: {}'.format(rf)

        if direction == 'f':
            seq  = fasta.seq
        else:
            seq = fasta.reverseComplement()

        pos = frame
        while pos < len(seq) - 2:
            codon = seq[pos:pos + 3]
            codon = codon.upper()
            # print('{}:{}:{}'.format(pos, codon, Fasta.codon2aa[codon]))
            trans.seq += Fasta.codon2aa[codon]
            pos += 3

        return trans

    def composition(fasta, uppercase=False):
        """-----------------------------------------------------------------------------------------
        Returns a dictionary with the composition of the sequence
        If uppercase is true, characters are converted to uppercase
        :return: dict, keys are letters in the sequence
        -----------------------------------------------------------------------------------------"""
        seq = fasta.seq
        if uppercase:
            seq = fasta.seq.upper()

        count = {}
        for ch in seq:
            if ch in count:
                count[ch] += 1
            else:
                count[ch] = 1

        return count

    def isACGT(fasta, threshold=0.8):
        """-----------------------------------------------------------------------------------------
        Return True if at least threshold fraction of characters in the sequence are ACGT

        :param threshold: float
        :return: Boolean
        -----------------------------------------------------------------------------------------"""
        total = fasta.length()
        if not total:
            return False

        comp = fasta.composition(uppercase=True)
        acgt = 0
        for base in 'ACGT':
            try:
                acgt += comp[base]
            except KeyError:
                # ignore missing bases
                continue

        if acgt / total > threshold:
            return True

        return False


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fasta = Fasta()
    fasta.id = 'sample1'
    fasta.doc = '20 each A,C,G,T'
    fasta.seq = 'A' * 20 + 'C' * 20 + 'G' * 20 + 't' * 20

    print(fasta.format(40))

    print('\nComposition')
    comp = fasta.composition()
    for ch in comp:
        print('\t{}\t{}'.format(ch, comp[ch]))
    comp = fasta.composition(uppercase=True)
    print('Uppercase')
    for ch in comp:
        print('\t{}\t{}'.format(ch, comp[ch]))

    # test isACGT
    print('\nACGT should be true')
    print(fasta.format(80))
    print('ACGT:{}'.format(fasta.isACGT()))

    print('\nACGT should be true')
    fasta.doc = 'No T'
    fasta.seq = 'A' * 20 + 'C' * 20 + 'G' * 20
    print(fasta.format(80))
    print('ACGT:{}'.format(fasta.isACGT()))

    print('\nACGT should be false')
    fasta.doc = 'ABC - 20 each'
    fasta.seq = 'A' * 20 + 'B' * 20 + 'C' * 20
    print(fasta.format(80))
    print('ACGT:{}'.format(fasta.isACGT()))

    # translation
    print('\nTranslation')
    fasta.id = 'sample1'
    fasta.doc = '20 each A,C,G,T'
    fasta.seq = 'A' * 20 + 'C' * 20 + 'G' * 20 + 't' * 20

    for direction in 'rf':
        for frame in range(3):
            trans = fasta.translate(frame=frame, direction=direction)
            print(trans.format(80))

exit(0)

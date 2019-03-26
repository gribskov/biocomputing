"""-------------------------------------------------------------------------------------------------
scaffold

15 October 2018     Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys


# import pysam


class Feature:
    """---------------------------------------------------------------------------------------------
    Sequence feature storage
    Not clear how to efficiently store this data so wrap it in a class to give independence from tne
    data representation

    code is a list of CIGAR codes
    countable is a list of codes we want to count
    refcountable is a list of CIGAR operations that advance the position in the reference sequence
    ---------------------------------------------------------------------------------------------"""
    code = ('match', 'insert', 'delete', 'skip', 'soft_clip', 'pad', 'equal', 'mismatch')
    countable = ('match', 'insert', 'delete', 'soft_clip', 'mismatch')
    refcountable = ('match', 'delete', 'skip', 'equal', 'mismatch')

    def __init__(self):
        self.alignment = None
        self.name = ''
        self.length = 0
        self.nreads = 0
        self.passed = 0
        self.feature = {}  # main data structure, hash of arrays of sequence position

    def get_alignment(self, bam, name=None):
        """-----------------------------------------------------------------------------------------
        sets the reference sequence and alignment bundle retrieved using pySAM.

        bam = pysam.AlignmentFile(bamfile, 'rb')    # pySAM AlignmentFile object
        for reference in bam.references:            # reference sequence name

        :param bam: pySAM AlignmentFile object
        :param name: reference,  sequence name (not tid)
        :return:
        -----------------------------------------------------------------------------------------"""
        if not name:
            return

        self.alignment = bam.fetch(name)
        self.name = name
        self.nreads = bam.count(name)
        self.length = bam.get_reference_length(name)

        return

    def from_cigar(self, read):
        """-----------------------------------------------------------------------------------------
        load features from CIGAR string
        cigar information is stored as tuples of (type, length)

        code    idx     meaning
        M       0       alignment match (can be a sequence match or mismatch)   BAM_CMATCH
        I       1       insertion to the reference                              BAM_CINS
        D       2       deletion from the reference                             BAM_CDEL
        N       3       skipped region from the reference                       BAM_CREF_SKIP
        S       4       soft clipping (clipped sequences present in SEQ)        BAM_CSOFT_CLIP
        H       5       hard clipping (clipped sequences NOT present in SEQ)    BAM_CHARD_CLIP
        P       6       padding (silent deletion from padded reference)         BAM_CPAD
        =       7       sequence match                                          BAM_CEQUAL
        X       8       sequence mismatch                                       BAM_CDIFF
        B       9       backwards skip                                          BAM_CBACK
        The B operator seems to have never been fully defined or implemented
        -----------------------------------------------------------------------------------------"""
        begin = read.reference_start

        # sys.stderr.write('from_cigar:{}\n'.format(read.query_name))
        if read.cigartuples:
            # None if not aligned
            # self.passed += 1
            for op in read.cigartuples:
                opcode = Feature.code[op[0]]
                size = op[1]
                # sys.stderr.write('     {} - {}     {}\n'.format(begin, begin + size, opcode))
                if opcode in Feature.countable:
                    # self.passed += 1
                    self.feature_add(opcode, begin, size, read)
                else:
                    sys.stderr.write(
                        'Feature::from_cigar - uncountable operation {}\n'.format(opcode))
                    sys.stderr.write('     {}:{} {} {}\n'.format(
                        self.name, read.query_name, begin, begin + size - 1))

                if opcode in Feature.refcountable:
                    begin += size

        else:
            sys.stderr.write('unmapped? {}\n'.format(read))

        return

    def feature_add(self, opcode, begin, size, read):
        """-----------------------------------------------------------------------------------------
        add counts to the feature count arrays
        -----------------------------------------------------------------------------------------"""
        feature = self.feature

        if not opcode in Feature.refcountable:
            # for features that do not advance the reference sequence posirion, just count 1
            size = 1

        try:
            self.passed += 1
            for i in range(begin, begin + size):
                # self.passed += 1
                feature[opcode][i] += 1
        except KeyError:
            feature[opcode] = [0 for _ in range(self.length)]
            # sys.stderr.write('Feature::feature-add - initializing {}\n'.format(opcode))
        except IndexError:
            sys.stderr.write(
                'Feature::feature-add - index error {}\t{}\t{}\n'.format(opcode, begin, size))
            sys.stderr.write('{}\n'.format(read))

        return

    def bases(self, opcode):
        """-----------------------------------------------------------------------------------------
        count of total bases in feature
        -----------------------------------------------------------------------------------------"""
        feature = self.feature

        bases = 0
        try:
            for i in range(len(feature[opcode])):
                bases += feature[opcode][i]

        except KeyError:
            bases = 0

        return bases


class Sam:
    """=============================================================================================
    Since Pysam is unreliable/hard to install, here is a simple parser for sam format
    ============================================================================================="""

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------
        it makes not sense to read the entire mapping file into meory as a default, so the class
        just allows youy to iterate over he file

        -----------------------------------------------------------------------------------------"""
        self.buffer = ''

        self.filename = filename
        self.fh = None
        if self.filename:
            self.opensam()

        self.header = []
        self.sequence = []
        self.program = []
        self.align = []

    def opensam(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        try:
            self.fh = open(self.filename, 'r')
        except IOError:
            sys.stderr.write('Sam::opensam - unable to open file ({})\n'.format(self.filename))
            exit(1)

        self.buffer = self.fh.readline().rstrip()

        return True

    def read_header(self):
        """-----------------------------------------------------------------------------------------
        Read and store header records from filehandle

        :return: int, number of header records
        -----------------------------------------------------------------------------------------"""
        header_n = 0
        while self.buffer.startswith('@HD'):
            self.header.append(self.buffer.rstrip())
            self.buffer = self.fh.readline()
            header_n += 1

        return header_n

    def read_sequence(self):
        """-----------------------------------------------------------------------------------------
        read and store sequence records from filehandle

        :return: int, number of sequence records
        -----------------------------------------------------------------------------------------"""
        seq_n = 0
        while self.buffer.startswith('@SQ'):
            self.sequence.append(self.buffer.rstrip())
            self.buffer = self.fh.readline()
            seq_n += 1

        return seq_n

    def read_program(self):
        """-----------------------------------------------------------------------------------------
        Read and store the program record

        :return:int, number of program records read
        -----------------------------------------------------------------------------------------"""
        pg_n = 0
        while self.buffer.startswith('@PG'):
            self.program.append(self.buffer.rstrip())
            self.buffer = self.fh.readline()
            pg_n += 1

        return pg_n

    def next(self):
        """-----------------------------------------------------------------------------------------
        Read the next alignment.  Assumes header, sequence, and program records have been read

        :return:
        -----------------------------------------------------------------------------------------"""
        field = self.buffer.split('\t')
        self.align = {
            'qname':field[0],
            'flag':int(field[1]),
            'rname':field[2],
            'pos': int(field[3]),
            'mapq': int(field[4]),
            'cigar':field[5],
            'rnext':field[6],
            'pnext':field[7],
            'tlen':int(field[8]),
            'seq':field[9],
            'qual':field[10]
        }
        self.buffer = self.fh.readline().rstrip()

        return field


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # old pysam based test

    # bamfile = sys.argv[1]
    # sys.stderr.write('bam file: {}\n'.format(bamfile))
    #
    # bam = pysam.AlignmentFile(bamfile, 'rb')

    # ref = {}
    # nref = 0
    # for i in range(len(bam.references)):
    #     ref[bam.references[i]] = bam.lengths[i]
    #     print( 'ref:{}     {}'.format(bam.references[i], bam.lengths[i]))
    #     nref += 1
    #     if nref > 99:
    #         break

    # sys.stderr.write('references: {}\n'.format(bam.nreferences))
    #
    # nbam = 0
    # nref = 0
    # for ref in bam.references:
    #     nref += 1
    #     if nref > 4:
    #         break
    #
    #     # p = 0
    #     # print(ref)
    #     # for pos in bam.pileup(ref):
    #     #     p += 1
    #     #     print(pos.n)
    #     #     if p > 5:
    #     # break
    #
    #     evidence = Feature()
    #     evidence.get_alignment(bam, ref)
    #     print('contig {}\n\t{} reads\n\t{} bases'.format(evidence.name, evidence.nreads,
    #                                                      evidence.length))
    #
    #     nread = 0
    #     for read in evidence.alignment:
    #         # print(read)
    #         nread += 1
    #         # sys.stderr.write('{}     {}\n'.format(nread, read.query_name))
    #         evidence.from_cigar(read)
    #
    #     print('contig {}\n\t{} reads\n\t{} bases'.format(evidence.name, evidence.nreads,
    #                                                      evidence.length))
    #     print('match: {}'.format(evidence.bases('match')))
    #     print('insert: {}'.format(evidence.bases('insert')))
    #     print('nref:', nref, '\tnread:', nread, '\tpassed', evidence.passed)

    samfile = sys.argv[1]
    sam = Sam(samfile)
    nheader = sam.read_header()
    print('{} header records'.format(nheader))

    nseq = sam.read_sequence()
    print('{} sequence records'.format(nseq))

    npg = sam.read_program()
    print('{} program records'.format(npg))

    mapped_n = 0
    link = {}
    while (sam.next()):
        mapped_n += 1
        rname = sam.align['rname']
        rnext = sam.align['rnext']

        if not rname in link:
            link[rname] = {}
            print('new r1 {}'.format(rname))

        if rnext in link[rname]:
            link[rname][rnext] += 1
        else:
            link[rname][rnext] = 1
            print('\tnew r2 {}'.format(rnext))

        if not mapped_n % 100000:
            print('{}\t{}'.format(mapped_n,sam.align['pos']))

        if not rname == 'NODE_1_length_118272_cov_3588.566145':
            break

    for r1 in link:
        for r2 in link[r1]:
            print('{}\t{}\t{}'.format(r1, r2, link[r1][r2]))

    for r1, r2 in

exit(0)

"""-------------------------------------------------------------------------------------------------
scaffold

15 October 2018     Michael Gribskov
-------------------------------------------------------------------------------------------------"""
import sys
import pysam


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


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    bamfile = sys.argv[1]
    sys.stderr.write('bam file: {}\n'.format(bamfile))

    bam = pysam.AlignmentFile(bamfile, 'rb')

    # ref = {}
    # nref = 0
    # for i in range(len(bam.references)):
    #     ref[bam.references[i]] = bam.lengths[i]
    #     print( 'ref:{}     {}'.format(bam.references[i], bam.lengths[i]))
    #     nref += 1
    #     if nref > 99:
    #         break

    sys.stderr.write('references: {}\n'.format(bam.nreferences))

    nbam = 0
    nref = 0
    for ref in bam.references:
        nref += 1
        if nref > 4:
            break

        # p = 0
        # print(ref)
        # for pos in bam.pileup(ref):
        #     p += 1
        #     print(pos.n)
        #     if p > 5:
        # break

        evidence = Feature()
        evidence.get_alignment(bam, ref)
        print('contig {}\n\t{} reads\n\t{} bases'.format(evidence.name, evidence.nreads,
                                                         evidence.length))

        nread = 0
        for read in evidence.alignment:
            # print(read)
            nread += 1
            # sys.stderr.write('{}     {}\n'.format(nread, read.query_name))
            evidence.from_cigar(read)

        print('contig {}\n\t{} reads\n\t{} bases'.format(evidence.name, evidence.nreads,
                                                         evidence.length))
        print('match: {}'.format(evidence.bases('match')))
        print('insert: {}'.format(evidence.bases('insert')))
        print('nref:', nref, '\tnread:', nread, '\tpassed', evidence.passed)

exit(0)

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
    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        self.alignment = None
        self.name = ''
        self.length = 0
        self.nreads = 0
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
        S       4       hard clipping (clipped sequences NOT present in SEQ)    BAM_CSOFT_CLIP
        P       6       padding (silent deletion from padded reference)         BAM_CPAD
        =       7       sequence match                                          BAM_CEQUAL
        X       8       sequence mismatch                                       BAM_CDIFF
        B       9       backwards skip                                          BAM_CBACK
        The B operator seems to have never been fully defined or implemented
        -----------------------------------------------------------------------------------------"""
        code = ('match', 'insert', 'delete', 'skip', 'soft_clip', 'pad', 'equal', 'mismatch')
        countable = ('match', 'insert', 'delete', 'soft_clip', 'mismatch')
        begin = read.reference_start

        sys.stderr.write('{}\n'.format(read.query_name))
        if read.cigartuples:
            # None if not aligned
            for op in read.cigartuples:
                opcode = code[op[0]]
                size = op[1]
                print('     {} - {}     {}     '.format(begin, begin + size, opcode))
                if opcode in countable:
                    self.feature_add(opcode, begin, size)
                else:
                    sys.stderr.write(
                        'Feature::from_cigar - uncountable operation {}\n'.format(opcode))
                    sys.stderr.write('     {}:{} {} {}\n'.format(
                        self.name, read.query_name, begin, begin + size - 1))

                begin += size
        # except TypeError:
        #     print('error:', read)

        return

    def feature_add(self, opcode, begin, size):
        """-----------------------------------------------------------------------------------------
        add counts to the feature count arrays
        -----------------------------------------------------------------------------------------"""
        pass

        return


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
    ncontig = 0
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
            block = read.get_blocks()
            evidence.from_cigar(read)
            # ct = read.cigartuples()
            # print(block,'/n',ct)
            if nread > 20:
                break

            # reference_start and reference are coordinates in ref

        print('ncontig:', ncontig)

exit(0)

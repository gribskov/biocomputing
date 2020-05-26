"""=================================================================================================
TRAPID produces lists of transcripts and their associated annotation terms, and annotation terms
with their asscociated transcripts.

To condense trinity transcript IDs, such as DN247_c1_g2_i3, to the gene level, it is easiest to
work with the second form: go_transcripts, interpro_transcripts, ko_transcripts etc.

Michael Gribskov     06 May 2020
================================================================================================="""
import sys


class Label:
    """---------------------------------------------------------------------------------------------
    Holds a single label and the list of sequences with the label
    ---------------------------------------------------------------------------------------------"""

    def __init__(self, line='', trim=0, sep='_'):
        self.n = 0
        self.id = ''
        self.description = ''
        self.evidence = ''
        self.member = {}

        if line:
            self.fromTrapidTerm(line, trim, sep)

    def addIds(self, idlist, trim=0, sep='_'):
        """----------------------------------------------------------------------------------------
         IDs to the label, if trim > 0, the ID string is split on sep, and trim number of tokens
         removed from the right end

         :param idlist: list, IDs of genes/transcipts
         :param trim: number of levels to remove from right side ID
         :param sep: string, separator character for ID string
         :return: int, number of IDs added
        -----------------------------------------------------------------------------------------"""
        for id in idlist:
            token = id.split(sep)
            for i in range(trim):
                token.pop()

            gene = sep.join(token)

            try:
                self.member[gene] += 1
            except KeyError:
                self.member[gene] = 1
                self.n += 1

        return self.n

    def fromTrapidTerm(self, line, trim, sep):
        """-----------------------------------------------------------------------------------------
        TRAPID produces lists of labels and the transcripts with the labels.  Each line is tab
        delimited with the following fields: counter, ID, evidence code, description, number of
        transcripts, transcriptlist

        the transcript list is a space delimeted string

        :param self:
        :param line: string, one line of file
        :return: string, label, ID field from line
        -----------------------------------------------------------------------------------------"""
        # field = line.split('\t')
        # counter = field.pop(0)
        # self.id = field.pop(0)
        # self.evidence = field.pop(0)
        # self.description = field.pop(0)
        # transcript_n = field.pop(0)
        (counter, self.id, self.evidence, self.description, transcript_n, transcripts) = line.split('\t')

        self.addIds(transcripts.split(' '), trim, sep)

        return self.id


class LabelSet():
    """---------------------------------------------------------------------------------------------
    Set of labels defined over a set of genes

    ---------------------------------------------------------------------------------------------"""

    def __init__(self, filename="", fh=None):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.term = {}
        self.gene = {}
        self.nterm = 0
        self.ngene = 0

        if fh:
            self.fh = fh
        elif filename:
            self.open(filename)

    def open(self, filename):
        """-----------------------------------------------------------------------------------------
        open a file for reading
        -----------------------------------------------------------------------------------------"""
        # print('opening:', filename)
        try:
            self.fh = open(filename, 'r')
        except (OSError, IOError):
            print('LabelSet::open - Unable to open file ({})'.format(filename))

        self.filename = filename
        return self.filename

    def fromTrapidTerm(self, trim=0, sep='_'):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        term_n = 0
        for line in self.fh:
            if line.startswith('#'):
                continue

            thislabel = Label(line, trim=trim, sep=sep)

            self.term[thislabel.id] = thislabel
            self.nterm += 1

            for gene in thislabel.member:
                try:
                    self.gene[gene] += 1
                except KeyError:
                    self.gene[gene] = 1
                    self.ngene += 1


        return self.nterm


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    ontoterms = LabelSet(filename=sys.argv[1])
    sys.stderr.write('\tTerm to transcript TRAPID file: {}\n'.format(ontoterms.filename))
    ontoterms.fromTrapidTerm(trim=1)
    print( '{} annotations for {} genes read from {}'.format(ontoterms.nterm, ontoterms.ngene,
                                                             sys.argv[1]))
    ref = open('ref.out', 'w')
    # write out as two columns, ID, term
    # sort by count of terms

    max = ontoterms.nterm / 3
    min = 5
    total = 0

    for term in sorted(ontoterms.term, key=lambda x: ontoterms.term[x].n, reverse=True ):
        if ontoterms.term[term].n > max:
            continue

        if ontoterms.term[term].n <min:
            break

        # print('{} {} {}'.format(term, ontoterms.term[term].n, ontoterms.term[term].description))
        for gene in ontoterms.term[term].member:
            ref.write('{} {}\n'.format(gene, term))
            total += 1

    print('{} annotations written'.format(total))
    ref.close()

# now the same thing for a selcted list
    selectout = open('select.out', 'w'
                  )
    selected_name = sys.argv[2]
    select = open(selected_name, 'r')
    selectlist = []
    for id in select:
        token = id.split('_')
        token.pop()
        gene = '_'.join(token)
        if gene not in selectlist:
            selectlist.append(gene)

    print('{} genes selected \n\n\n'.format(len(selectlist)))

    total = 0
    for term in sorted(ontoterms.term, key=lambda x: ontoterms.term[x].n, reverse=True):
        if ontoterms.term[term].n > max:
            continue

        if ontoterms.term[term].n < min:
            break

        # print('{} {} {}'.format(term, ontoterms.term[term].n, ontoterms.term[term].description))
        for gene in ontoterms.term[term].member:
            if gene in selectlist:
                selectout.write('{} {}\n'.format(gene, term))
                total += 1

    print('{} selected annotations written'.format(total))
    selectout.close()

    exit(0)

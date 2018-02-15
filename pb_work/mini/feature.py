class Feature:
    """
    The Feature class is for describing ranges in a linear coordinate system, such as are described
    in GFF3 or GTF files.
    """

    def __init__(self):

        self.references = []  # references for coordinates, e.g. chromosomes
        self.features = []  # individual features

        return None

    def readGFF3(self, filename, feature_type=''):
        """
        Read a set of features from a GFF3 file. Example:
        Pt	ensembl	protein_coding_gene	15938	20068	.	-	.	ID=ATCG00170;biotype=protein_coding;description=DNA-directed RNA polymerase family protein [Source:TAIR_LOCUS%3BAcc:ATCG00170];external_name=RPOC2;logic_name=tair

        Columns:
            0: seqid (reference sequence)
            1: source
            2: feature
            3: begin
            4: end
            5: score (float, meaning varies)
            6: strand
            7: phase (i.e., reading frame)
            8: attributes

        :param filename: name of filename to read
        :param feature_type: name of feature to read, default = '', i.e., all
        :return: number of features
        """
        try:
            gffin = open(filename, 'r')
        except Exception as err:
            print('Features::readGFF - unable to read GFF file ({})'.format(filename))
            print(err)
            exit(1)

        columns = []
        nfeature = 0
        for line in gffin:
            if line.isspace() or line.startswith('#'):
                continue

            if not feature_type in line:
                # select only matching features
                continue

            col = line.split()
            feature = {'seqid': col[0],
                       'source': col[1],
                       'feature': col[2],
                       'begin': col[3],
                       'end': col[4],
                       'score': col[5],
                       'strand': col[6],
                       'phase': col[7],
                       }
            if feature['strand'] == '.' or feature['strand'] == '?':
                feature['strand'] = None
            if feature['score'] == '.':
                feature['score'] = None

            attrs = col[8].split(';')
            attribute = {}
            for att in attrs:
                left, right = att.split('=')
                attribute[left] = right
            feature['attribute'] = attribute
            self.features.append(feature)
            nfeature += 1

        return nfeature

        # end of readGFF3


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gff = Feature()
    n = gff.readGFF3('at_1000k.gff3', feature_type='CDS')
    print('{} features read'.format(n))

    exit(0)

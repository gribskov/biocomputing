class Feature:
    """=============================================================================================
    The Feature class is for describing ranges in a linear coordinate system, such as are described
    in GFF3 or GTF files.

    Synopsis
        # Read a GFF file
        gff = Feature()
        gff.label = 'gff'
        gff_file = 'at_1000k.gff3'
        n = gff.readGFF3(gff_file, feature_type='gene')
        
        $ Read a blast tabular output with evalue cutoff of 1e-40
        blast = Feature()
        blast.label = 'blast'
        blast_file = 'ch4.blastn'
        n = blast.readBlastTabular(blast_file, evalue_cutoff=1e-40)
    
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Feature constructor:
            Feature is a collection of individual features.

            label: is a string identifying the source of the features
            references: list, the sequences on which features may be located
            features: list of dictionary - fields vary with the datatype.
                see readGFF3 and readBlastTabular
        -----------------------------------------------------------------------------------------"""

        self.references = []  # reference sequences for coordinates, e.g. chromosomes
        self.features = []  # individual features
        self.label = ''

    def readGFF3(self, filename, feature_type=''):
        """-----------------------------------------------------------------------------------------
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
        -----------------------------------------------------------------------------------------"""
        gffin = None
        try:
            gffin = open(filename, 'r')
        except Exception as err:
            print('Features::readGFF - unable to read GFF file ({})'.format(filename))
            print(err)
            exit(1)

        nfeature = 0
        for line in gffin:
            if line.isspace() or line.startswith('#'):
                continue

            col = line.split()
            if feature_type not in col[2]:
                # select only matching features
                continue

            feature = {'seqid': col[0],
                       'source': col[1],
                       'feature': col[2],
                       'begin': int(col[3]),
                       'end': int(col[4]),
                       'score': col[5],
                       'strand': col[6],
                       'phase': col[7],
                       }
            if feature['strand'] == '.' or feature['strand'] == '?':
                feature['strand'] = None

            if feature['score'] == '.':
                feature['score'] = 0.0
            else:
                # if score is not '.' it should be a float
                feature['score'] = float(feature['score'])

            if not feature['seqid'] in self.references:
                # keep a list of the reference sequences
                self.references.append(feature['seqid'])

            # the variable attributes in column 9
            attrs = col[8].split(';')
            for att in attrs:
                left, right = att.split('=')
                feature[left] = right

            if 'Parent' in feature:
                # for pseudogene child features
                feature['ID'] = feature['Parent']

            self.features.append(feature)
            nfeature += 1

        return nfeature

        # end of readGFF3

    def readBlastTabular(self, filename, evalue_cutoff=1e-40):
        """-----------------------------------------------------------------------------------------
        Read blast -m 8 tabular format. 12 columns:
        0	Query	The query sequence id
        1	Subject	The matching subject sequence id
        2	% id
        3	alignment length
        4	mismatches
        5	gap openings
        6	q.start
        7	q.end
        8	s.start
        9	s.end
        10	e-value
        11	bit score

        :param filename: file to read
        :param: evalue_cutoff - maximum E-value
        :return: nfeature
        -----------------------------------------------------------------------------------------"""
        blastin = None
        try:
            blastin = open(filename, 'r')
        except Exception as err:
            print('Features::readBlastTabular - unable to read Blast file ({})'.format(filename))
            print(err)
            exit(1)

        nfeature = 0
        for line in blastin:
            if line.isspace() or line.startswith('#'):
                continue

            field = line.rstrip().split()
            if float(field[10]) > evalue_cutoff:
                continue

            # make sure begin is less than end
            begin = int(field[8])
            end = int(field[9])
            if end < begin:
                begin, end = end, begin

            self.features.append({'seqid': field[1], 'begin': begin, 'end': end, 'ID': field[0]})
            nfeature += 1

        return nfeature

        # end of readBlastTabular

    def sortByPos(self):
        """-----------------------------------------------------------------------------------------
        Sorts the features by seqid and beginning position.  Sort is done in place.

        :return: None
        -----------------------------------------------------------------------------------------"""
        self.features.sort(key=lambda f: (f['seqid'], f['begin']))

        return None

    def ranges(self):
        """-----------------------------------------------------------------------------------------
        Merges overlapping features, feature order may change due to sorting by position. because
        IDs may overlap, for instant, seq1 and seq11, the id list is terminated with ;
        The final ; must be removed when the range closes.
            current['ID'] = current['ID'].rstrip(';')

        :return: feature object with ranges merged
        -----------------------------------------------------------------------------------------"""
        self.sortByPos()
        range = Feature()
        range.label = self.label

        # The first feature always starts a new range
        range.features.append(
            {'seqid': self.features[0]['seqid'],
             'feature': 'range',
             'begin': self.features[0]['begin'],
             'end': self.features[0]['end'],
             'ID': self.features[0]['ID'] + ';'
             })

        # for convenience, make a pointer to most recent feature
        current = range.features[-1]

        for feature in self.features:

            if feature['seqid'] != current['seqid'] or feature['begin'] > current['end']:
                # new sequence always requires new range
                # otherwise, start a new range if begin <= end
                current['ID'] = current['ID'].rstrip(';')
                range.features.append(
                    {'seqid': feature['seqid'],
                     'feature': 'range',
                     'begin': feature['begin'],
                     'end': feature['end'],
                     'ID': feature['ID'] + ';'
                     })
                current = range.features[-1]

            else:
                # overlap: continue previous range i - update end and id
                current['end'] = max(current['end'], feature['end'])
                if feature['ID'] + ';' not in current['ID']:
                    # add id if not present if id list
                    current['ID'] += '{};'.format(feature['ID'])

        current['ID'] = current['ID'].rstrip(';')

        return range

        # end of ranges


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # Read the GFF file
    gff = Feature()
    gff_file = 'at_1000k.gff3'
    n = gff.readGFF3(gff_file, feature_type='gene')
    gff.label = 'gff'
    print('Read gff annotation')
    print('     {} features read from {}'.format(n, gff_file))

    print('     reference sequences found in {}:'.format(gff_file), end='')
    for id in sorted(gff.references):
        print(' {}'.format(id), end='')
    print()

    # find the overlapping GFF ranges and write to output file
    out = open('gffranges.test.txt', 'w')
    n = 0
    gffranges = gff.ranges()
    for r in gffranges.features:
        out.write(
            '{} {} {} {} {}\n'.format(r['feature'], r['seqid'], r['begin'], r['end'], r['ID']))
        n += 1
    print('     {} ranges found in {}'.format(n, gff_file))
    out.close()

    # Read the blast tabular output
    print('\nRead blast tabular')
    blast_file = 'ch4.blastn'
    blast = Feature()
    blast.label = 'blast'
    n = blast.readBlastTabular(blast_file, evalue_cutoff=1e-40)
    print('    {} sequences read from {}'.format(n, blast_file))

    # Find the overlapping blast ranges and write to output file
    n = 0
    blastranges = blast.ranges()
    out = open('blastranges.test.txt', 'w')
    for r in blastranges.features:
        out.write(
            '{} {} {} {} {}\n'.format(r['feature'], r['seqid'], r['begin'], r['end'], r['ID']))
        n += 1

    print('    {} ranges found in {}'.format(n, blast_file))
    out.close()

    exit(0)

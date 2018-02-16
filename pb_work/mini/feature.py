class Feature:
    """=============================================================================================
    The Feature class is for describing ranges in a linear coordinate system, such as are described
    in GFF3 or GTF files.
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Feature constructor
        -----------------------------------------------------------------------------------------"""

        self.references = []  # references for coordinates, e.g. chromosomes
        self.features = []  # individual features

    def __iter__(self):
        """-----------------------------------------------------------------------------------------
        Iterator is a generator
        :return: generator
        -----------------------------------------------------------------------------------------"""
        return self.next()

    def readGFF3(self, filename, feature_type=''):
        """-----------------------------------------------------------------------------------------
        Read a set of features from a GFF3 file. Example:
        Pt	ensembl	protein_coding_gene	15938	20068	.	-	.	ID=ATCG00170;biotype=protein_coding;description=DNA-directed RNA polymerase family protein [Source:TAIR_LOCUS%3BAcc:ATCG00170];external_name=RPOC2;logic_name=tair
        TODO: bug in selection, feature_type matches in feature including the key
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

            if feature_type not in line:
                # select only matching features
                continue

            col = line.split()
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

            nfeature += 1

        return nfeature

        # end of readGFF3

    def next(self):
        """-----------------------------------------------------------------------------------------
        generator for iterating over features
        -----------------------------------------------------------------------------------------"""
        for feature in self.features:
            yield feature
        raise StopIteration

    def sortByPos(self):
        """-----------------------------------------------------------------------------------------
        Sorts the features by seqid and beginning position.  Sort is done in place.

        :return: None
        -----------------------------------------------------------------------------------------"""
        self.features.sort(key=lambda f: (f['seqid'], f['begin']))

        return None

    def ranges(self):
        """-----------------------------------------------------------------------------------------
        Merges overlapping ranges, feature order may change due to sorting by position

        :return: feature object with ranges merged
        -----------------------------------------------------------------------------------------"""
        self.sortByPos()
        range = Feature()

        # start a range feature
        begin = 0
        end = 0
        seqid = ''
        idlist = []

        for feature in self.features:
            if feature['seqid'] == seqid:
                if feature['begin'] <= end:
                    # continue previous range
                    end = feature['end']
                    begin = feature['begin']
                    if 'ID' in feature:
                        if feature['ID'] not in idlist:
                            idlist.append['ID'])
                else:
                # save current
                range.features.append({'seqid': seqid, 'begin': begin, 'end': end, 'ID': ';'.join(idlist)})

                # start new range
                begin = feature['begin']
                end = feature['end']
                idlist = []
                if 'ID' in feature:
                    idlist.append(feature['ID'])

            else:
                # new sequence, must start new range

                # save current
                    range.features.append({'seqid': seqid, 'begin': begin, 'end': end, 'ID': ';'.join(idlist)})

                # new range
                begin = feature['begin']
                end = feature['end']
                idlist = []
                if 'ID' in feature:
                    idlist.append(feature['ID'])

        return range



# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gff = Feature()
    n = gff.readGFF3('at_1000k.gff3', feature_type='gene')
    print('{} features read'.format(n))

    print('reference sequences:')
    for id in gff.references:
        print('    {}'.format(id))

    print('\niteration with generator')
    n = 0
    for f in gff.next():
        print('1:{} {} {} {}'.format(f['seqid'], f['begin'], f['end'], f['feature']))
        n += 1
        if n > 5:
            break

    print('\niteration with iterator, after sort by pos')
    gff.sortByPos()
    n = 0
    for f in gff:
        print('2:{} {} {} {}'.format(f['seqid'], f['begin'], f['end'], f['feature']))
        print('    {}'.format(f['attribute']))
        n += 1
        if n > 5:
            break

    exit(0)

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
        self.lable = ''

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
        4	mistmatches
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
        range.label = self.label

        # start a range feature
        begin = self.features[0]['begin']
        end = self.features[0]['end']
        seqid = self.features[0]['seqid']
        idlist = []

        for feature in self.features:
            if feature['seqid'] == seqid:
                if feature['begin'] <= end:
                    # continue previous range
                    end = max(end, feature['end'])

                else:
                    # save current
                    range.features.append(
                        {'seqid': seqid, 'feature': 'range', 'begin': begin, 'end': end,
                         'ID': ';'.join(idlist)})

                    # start new range
                    begin = feature['begin']
                    end = feature['end']
                    seqid = feature['seqid']
                    idlist = []

                if 'ID' in feature:
                    if feature['ID'] not in idlist:
                        idlist.append(feature['ID'])

            else:
                # new sequence, must start new range

                # save current
                range.features.append(
                    {'seqid': seqid, 'feature': 'range', 'begin': begin, 'end': end,
                     'ID': ';'.join(idlist)})

                # new range
                begin = feature['begin']
                end = feature['end']
                seqid = feature['seqid']
                idlist = []
                if 'ID' in feature:
                    if feature['ID'] not in idlist:
                        idlist.append(feature['ID'])

        return range


    def compareRange(self, other):
        """---------------------------------------------------------------------------------------------
        must define seqid, begin, end
    
        :param self: first range
        :param other: second range
        :return: overlap
        ---------------------------------------------------------------------------------------------"""
        print('first:', self.features[self.ptr])
        print('second:', other.features[other.ptr])

        first = self.features[self.ptr]
        second = other.features[other.ptr]

        if first.seqid < second.seqid or first.end < second.begin:
            return 'left'
        elif first.seqid > second.seqid or first.begin > second.end:
            return 'right'
        else:
            return 'overlap'


        return

def compareFeature(f1,f2):
    """---------------------------------------------------------------------------------------------
    returned feature has begin <= to other
    :param f1:
    :param f2:
    :return:
    ---------------------------------------------------------------------------------------------"""
    if f1['seqid'] < f2['seqid']:
        return f1
    elif f1['seqid'] > f2['seqid']:
        return f2
    
    # seqid are equal

    if f1['begin'] <= f2['begin']:
        return f1
    elif f1['begin'] > f2['begin']:
        return f2

# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    LIMIT = 1000

    gff = Feature()
    n = gff.readGFF3('at_1000k.gff3', feature_type='gene')
    gff.label = 'gff'
    print('{} features read'.format(n))

    print('reference sequences:')
    for id in gff.references:
        print('    {}'.format(id))

    # print('\niteration with generator')
    # n = 0
    # for f in gff.next():
    #     print('1:{} {} {} {}'.format(f['seqid'], f['begin'], f['end'], f['feature']))
    #     n += 1
    #     if n > LIMIT:
    #         break

    print('\niteration with iterator, after sort by pos')
    gff.sortByPos()
    n = 0
    for f in gff:
        print('range:{} {} {} {}'.format(f['seqid'], f['begin'], f['end'], f['ID']))
    n += 1
    # if n > LIMIT:
    #     break

    print('\n    ranges')
    out = open('gffranges.mrg.txt', 'w')
    n = 0
    gffranges = gff.ranges()
    for r in gffranges.features:
        out.write(
            '{} {} {} {} {}\n'.format(r['feature'], r['seqid'], r['begin'], r['end'], r['ID']))
    print('    range:{} {} {} {}'.format(r['seqid'], r['begin'], r['end'], r['ID']))
    n += 1
    # if n > LIMIT:
    #     break

    out.close()

    print('\nread blast tabular')
    blast = Feature()
    blast.label = 'blast'
    n = blast.readBlastTabular('ch4.blastn', evalue_cutoff=1e-40)
    print('    {} sdequences read'.format(n))
    print('\n    ranges')
    n = 0
    blastranges = blast.ranges()
    out = open('blastranges.mrg.txt', 'w')
    for r in blastranges.features:
        out.write(
            '{} {} {} {} {}\n'.format(r['feature'], r['seqid'], r['begin'], r['end'], r['ID']))
    print('    range:{} {} {} {}'.format(r['seqid'], r['begin'], r['end'], r['ID']))
    n += 1
    # if n > LIMIT:
    #     break

    out.close()

    joint = Feature()
    gffranges.ptr = 0
    blastranges.ptr = 0
    overlap = gffranges.compareRange(blastranges)
    if overlap == 'left':
        joint.features.append(gffranges.features[gffranges.ptr])
        gffranges.ptr += 1
    elif overlap == 'right'
        joint.features.append(blastranges.features[blastranges.ptr])
        blastranges.ptr += 1
    else:
        # TODO process overlap
        pass

    exit(0)

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
        Feature constructor
        -----------------------------------------------------------------------------------------"""

        self.references = []  # references for coordinates, e.g. chromosomes
        self.features = []  # individual features
        self.label = ''

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
        Merges overlapping features, feature order may change due to sorting by position. because
        IDs may overalp, for instant, seq1 and se11, the id list is terminated with ;
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
             'ID': self.features[0]['ID']+';'
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
                     'ID': feature['ID']+';'
                     })
                current = range.features[-1]

            else:
                # overlap: continue previous range i - update end and id
                current['end'] = max(current['end'], feature['end'])
                if feature['ID']+';' not in current['ID']:
                    # add id if not present if id list
                    current['ID'] += '{};'.format(feature['ID'])

        current['ID'] = current['ID'].rstrip(';')

        return range

        # end of ranges

    def compareRange(self, other):
        """-----------------------------------------------------------------------------------------
        The feature objects, self and other, must define a label attribute.
        The ranges, i.e., the features within each object, must define seqid, begin, end

        The returned Feature object defines features with: label, seqid, begin, end
        the label is one of the input labels, or 'both' for segments found in both
    
        :param self: Feature object, the first set of ranges
        :param other: Feature object, the second set of ranges
        :return: Feature object
        -----------------------------------------------------------------------------------------"""
        joint = Feature()
        # add pointers to features list
        self.ptr = 0
        other.ptr = 0

        a = nextRange(self, other)
        while True:

            b = nextRange(self, other)
            if not b:
                break

            # nextRange ensures that a will always be less than b,either with smaller seqid, if a is
            # on an earlier sequence (chromosome), or a smaller begin if on the same chromosome.

            # if a ends before b begins, accept the whole a range
            if a['seqid'] < b['seqid'] or a['end'] < b['begin']:
                joint.features.append({'label': a['label'], 'seqid': a['seqid'],
                                       'begin': a['begin'], 'end': a['end']})
                a = b
                continue

            # a and b overlap, three possible segments
            #   a segment in a only (missing if a['begin'] == b['begin'])
            #   a segment in both (always present because a and b overlap)
            #   an overhang at the end (maybe missing if a['end'] == b['end']
            #       this segment may overlap the next range so it is saved as a for the next cycle

            begin = a['begin']
            if a['begin'] < b['begin']:
                # a only segment
                joint.features.append({'label': a['label'], 'seqid': a['seqid'],
                                       'begin': a['begin'], 'end': b['begin'] - 1})

            begin = b['begin']
            if a['end'] >= b['end']:
                # both segment, b ends first
                joint.features.append({'label': 'both', 'seqid': a['seqid'],
                                       'begin': begin, 'end': b['end']})
                a['begin'] = b['end'] + 1
                # it is possible possible that a['end'] is < a['begin'].   this is ok because in the
                # next cycle the segment will not overlap
                continue

            if a['end'] < b['end']:
                # both segment, a ends first
                joint.features.append({'label': 'both', 'seqid': a['seqid'],
                                       'begin': begin, 'end': a['end']})
                b['begin'] = a['end'] + 1
                a = b
                continue

        return joint

        # end of compareRange


# non-object functions -----------------------------------------------------------------------------

def nextRange(r1, r2):
    """---------------------------------------------------------------------------------------------
    return the next feature in ranges r1 and r2
    :param r1: Feature object
    :param r2: Feature object
    :return: Feature.features element
    ---------------------------------------------------------------------------------------------"""
    # cases where one or both ranges have been completely processed
    if r1.ptr >= len(r1.features) and r2.ptr >= len(r2.features):
        # both ranges are completed
        return None

    elif r1.ptr >= len(r1.features):
        # range 1 only is completed
        f2 = r2.features[r2.ptr]
        f2['label'] = r2.label
        r2.ptr += 1
        return f2

    elif r2.ptr >= len(r2.features):
        # range 2 only is completed
        f1 = r1.features[r1.ptr]
        f1['label'] = r1.label
        r1.ptr += 1
        return f1

    # general case where both ranges have available features    
    f1 = r1.features[r1.ptr]
    f2 = r2.features[r2.ptr]
    f1['label'] = r1.label
    f2['label'] = r2.label

    if f1['seqid'] < f2['seqid']:
        r1.ptr += 1
        return f1

    elif f1['seqid'] > f2['seqid']:
        r2.ptr += 1
        return f2

    elif f1['begin'] <= f2['begin']:
        r1.ptr += 1
        return f1

    else:
        r2.ptr += 1
        return f2

    # end of nextRange


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
    out = open('gffranges.mrg.txt', 'w')
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
    out = open('blastranges.mrg.txt', 'w')
    for r in blastranges.features:
        out.write(
            '{} {} {} {} {}\n'.format(r['feature'], r['seqid'], r['begin'], r['end'], r['ID']))
        n += 1

    print('    {} ranges found in {}'.format(n, blast_file))
    out.close()

    # Find the overlaps between the GFF and Blast ranges
    joint = gffranges.compareRange(blastranges)

    # write overlapping segments to a file
    filename = 'overlap.mrg.txt'
    try:
        out = open(filename, 'w')
    except Exception as err:
        print('Unable to open output file ({})'.format(filename))
        print(err)
        exit(1)

    # print a report of the overlap between the GFF and Blast ranges
    count = {'total': 0}
    for overlap in joint.features:
        out.write('{}\t{}\t{}\t{}\n'.
                  format(overlap['label'], overlap['seqid'], overlap['begin'], overlap['end']))
        try:
            count[overlap['label']] += overlap['end'] - overlap['begin'] + 1
        except KeyError:
            count[overlap['label']] = overlap['end'] - overlap['begin'] + 1
        count['total'] += overlap['end'] - overlap['begin'] + 1

    print('\nCounts by label')
    for type in count:
        print('{:5}{:10d}'.format(type, count[type]))

    print('\n{:5}{:>10}{:>8}{:>8}'.format('label', 'count', '%total', '%both'))
    for type in ['gff', 'blast']:
        print('{:5}{:10d}{:8.2f}{:8.2f}'.
              format(type, count[type],
                     100.0 * count[type] / count['total'], 100.0 * count['both'] / count[type],
                     ))

    out.close()

    exit(0)

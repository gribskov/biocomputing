"""=================================================================================================
compare data described as linear ranges on a line, e.g., genes on a chromosome

This code was developed for the practical bioinformatics class, and was formerly called
genome_ranges. I changed the name to not cause confusion with the R package Genomic_Ranges

21 March 2025   Michael Gribskov
================================================================================================="""


class Range:
    """---------------------------------------------------------------------------------------------
    The Range class is for describing ranges on a linear coordinate system, such as are
    described in GFF3/GTF or blast result files.

    Synopsis
        # Read a GFF file

        gff = Range('at_1000k.gff3')
        n = gff.read_gff(select=['gene','transcript'])

        $ Read a blast tabular output with evalue cutoff of 1e-40

        blast = Range('ch4.blastn')
        n = blast.readBlastTabular(evalue_cutoff=1e-40)
        
    See Testing section at bottom file for usage examples
    ---------------------------------------------------------------------------------------------"""

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------
        Range constructor:
        
        filename, filehandle: path to data file, and open filehandle
        features: list of dictionary.  keys:
            sequence: chromosome or scaffold name
            begin: beginning position of sequence region, zero based
            end: ending position of sequence region
            name: semicolon delimited strings of entities contributing to the range, i.e., the
                original transcripts or blast hits

        :param filename: string, file with input data
        -----------------------------------------------------------------------------------------"""

        self.filename = filename  # name of file
        self.file = None  # filehandle
        self.features = []  # list of selected features

        if filename:
            self.open_file('r')

    def open_file(self, mode):
        """-----------------------------------------------------------------------------------------
        Opens the file in self.filename with error checking and sets the filehandle, self.file
        Exit with status=1 if the file can't be opened

        param mode: string, mode for opening file, usually 'r' or 'w'
        return: filehandle of opened file
        -----------------------------------------------------------------------------------------"""
        try:
            self.file = open(self.filename, mode)
        except OSError:
            print(f'Error opening file ({self.filename})')
            # since you cannot continue, exit with status = 1
            exit(1)

        return self.file

    def read_gff(self, select=None):
        """-----------------------------------------------------------------------------------------
        Populate the features with a list dictionaries from a GFF file. The dictionary has the keys
        sequence, name, begin, end, where sequence is the ID of the genomic sequences and name is
        the feature ID
        Rows are selected if the feature (column 3)  is in the parameter select.

        example GFF format:
        scaffold9   AUGUSTUS    gene        135854  136842  .   +   .   ID=jg677
        scaffold9   AUGUSTUS    transcript  135854  136842  .   +   .   ID=jg677.t1;Parent=jg677
        scaffold9   AUGUSTUS    start_codon 135854  135856  .   +   0   ID=jg677.t1.start;Parent=jg677.t1
        scaffold9   AUGUSTUS    CDS         135854  135929  1   +   0   ID=jg677.t1.e0;Parent=jg677.t1
        sequence    source      feature     begin   end  score strand frame attributes

        param select: list of strings corresponding to GFF fields
        return: int, number of features read
        -----------------------------------------------------------------------------------------"""
        if not select:
            select = ['transcript']

        n_feature = 0
        for line in self.file:
            field = line.split('\t')
            if len(field) < 9:
                # skip records with less than 9 fiels (not data records)
                continue

            # all features
            gff = {'sequence': field[0], 'source': field[1], 'feature': field[2],
                   'begin':    int(field[3]), 'end': int(field[4]), 'strand': field[6],
                   'frame':    field[7], 'attr': field[8]}

            # the attributes are a semicolon delimited list of tag value pairs; the feature ID is in the attributes
            attrs = gff['attr'].split(';')
            attribute = {}
            for a in attrs:
                tag, value = a.split('=')
                attribute[tag] = value
            if not attribute['ID']:
                print(attribute)

            if gff['feature'] in select:
                self.features.append({'sequence': gff['sequence'], 'name': attribute['ID'],
                                      'begin':    gff['begin'], 'end': gff['end']})
                n_feature += 1

        return n_feature

    def read_blast_tabular(self, e_threshold=1e-50):
        """-------------------------------------------------------------------------------------------------------------
        Returns a list of features from a BlastX results, such as:
            - scaffold number, start position, end position, subject/geneID
        -------------------------------------------------------------------------------------------------------------"""
        self.features = []  # remove any existing features

        for line in self.file:
            if not line:
                # end of file
                break

            if line.isspace():
                # blank line
                continue

            field = line.rstrip().split('\t')
            if float(field[10]) < 0:
                # skip damaged records with impossible e-values
                continue

            begin = int(field[6])
            end = int(field[7])
            if begin > end:
                # make sure begin is less than end
                begin, end = end, begin

            subject = field[1].replace('UniRef50_', '')
            ul_pos = field[0].index('_')
            query = field[0][:ul_pos]

            if float(field[10]) <= e_threshold:
                self.features.append({'sequence': query, 'begin': begin, 'end': end, 'name': subject})

        return len(self.features)

    def sort_by_pos(self):
        """-----------------------------------------------------------------------------------------
        Sorts the features by sequence and beginning position.  Sort is done in place.

        usage
            self.sort_by_pos()

        :return: None
        -----------------------------------------------------------------------------------------"""
        self.features.sort(key=lambda f: (f['sequence'], f['begin']))

    def overlap(self):
        """-----------------------------------------------------------------------------------------
        Merges overlapping features, feature order may change due to sorting by position. because
        IDs may overlap, for instant, seq1 and seq11, the id list is terminated with ;
        The final ; must be removed when the range closes.
            current['name'] = current['name'].rstrip(';')

        :return: new feature object with ranges merged
        -----------------------------------------------------------------------------------------"""
        self.sort_by_pos()
        lrange = Range()  # output overlapped range object

        # The first feature always starts a new range
        feature = {'sequence': self.features[0]['sequence'],
                   'name':     self.features[0]['name'],
                   'begin':    self.features[0]['begin'],
                   'end':      self.features[0]['end']}

        lrange.features.append(feature)

        # for convenience, make a pointer to most recent feature
        current = lrange.features[-1]

        for feature in self.features[1:]:

            if feature['sequence'] != current['sequence'] or feature['begin'] > current['end'] + 1:
                # new sequence always requires new range
                # otherwise, start a new range if begin > end + 1
                # a range that ends at n and a feature that begins at n+1 are considered overlapped
                lrange.features.append({'sequence': feature['sequence'],
                                       'name':     feature['name'],
                                       'begin':    feature['begin'],
                                       'end':      feature['end']})
                current = lrange.features[-1]

            else:
                # overlap: continue previous range i - update end and id
                current['end'] = max(current['end'], feature['end'])
                current['name'] = current['name'] + f';{feature["name"]}'

        return lrange

    def make_list(self, feature, maxnames=3):
        """-----------------------------------------------------------------------------------------
        sort the name list and cut off at maxnames explicitly displayed names.  Extra names are
        shown with the phrase ' + n more'

        :param maxnames: int, maximum number of explicit names
        :param feature: dict, keys = ['sequence', 'name', 'begin', 'end' ]
        :return: string, formatted string for printing
        -----------------------------------------------------------------------------------------"""
        names = feature['name'].split(';')
        if len(names) > maxnames:
            more = len(names) - maxnames
            namelist = ';'.join(sorted(names)[:maxnames]) + f' + {more} more'
        else:
            namelist = ';'.join(sorted(names))

        return namelist


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # Read the GFF file
    gff_file = 'Hp.augustus.hints.gff'
    gff = Range(filename=gff_file)
    n = gff.read_gff(select=['transcript'])
    gff.label = 'gff'
    print('Read gff annotation')
    print(f'\t{n} features read from {gff_file}')

    # find the overlapping GFF ranges and write to output file
    gffranges = gff.overlap()
    gffranges.filename = 'gffranges.test.txt'
    gffranges.open_file('w')

    n = 0
    for r in gffranges.features:
        namelist = gffranges.make_list(r, 3)
        gffranges.file.write(f'range {r["sequence"]} {r["begin"]} {r["end"]} {namelist}\n')
        n += 1
    print(f'\t{n} ranges found in {gff_file}')
    gffranges.file.close()

    # Read the blast tabular output
    print('\nRead blast tabular')
    blast_file = 'genome.small_uniref50.dmndblastx'
    blast = Range(filename=blast_file)
    blast.label = 'blast'
    nblast = blast.read_blast_tabular(e_threshold=1e-50)
    print('    {} sequences read from {}'.format(nblast, blast_file))

    # Find the overlapping blast ranges and write to output file
    blast.sort_by_pos()
    blast_range = blast.overlap()
    blast_range.filename = 'blastranges.test.txt'
    blast_range.open_file('w')

    n = 0
    for f in blast_range.features:
        namelist = blast_range.make_list(f, 3)
        blast_range.file.write(f"range {f['sequence']} {f['begin']} {f['end']} {namelist}\n")
        n += 1

    print('    {} ranges found in {}'.format(n, blast_file))
    blast_range.file.close()

    exit(0)

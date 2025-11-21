"""#################################################################################################
Another GFF/GTF (GXF) parser with the individual records and collections of records split into two
classes. This should make it easier reuse.

Michael Gribskov     21 November 2025
#################################################################################################"""
version = '2.0.0'


class GxfRecord:
    """=============================================================================================
    The GFF standard defines columns 0-7, column 8 is used for any attributes the project defines
    For details see https://useast.ensembl.org/info/website/upload/gff.html or
    https://agat.readthedocs.io/en/latest/gxf.html

    Column 0: "seqid"
      The ID of the landmark used to establish the coordinate system for the current feature,
      typically the ID of the chromosome or scaffold. IDs may contain any characters, but must
      escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not
      contain unescaped whitespace and must not begin with an unescaped ">".
    Column 1: "source"
      Free text describing the origin of the feature annotation. Typically this is the name of
      program, such as "Genescan" or a database name, such as "Genbank."
    Column 2: "type"
      The feature type. This is constrained to be either a term from the Sequence Ontology. Common
      features include gene, transcript, mRNA, 5'UTR, 3'UTR, exon, CDS, start_codon, stop_codon, etc
    Columns 3 & 4: "start" and "end"
      The start and end coordinates of the feature in 1-based integer coordinates, relative to the
      sequence in column 1. Start is always less than or equal to end. For features that cross the
      origin of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes),
      end end + the length of sequence, i.e., the end will be past the end of the sequence.
      For zero-length features, such as insertion sites, start equals end and the implied site is to
      the right of the indicated base in the direction of the landmark.
    Column 5: "score"
      The score of the feature, a floating point number. As in earlier versions of the format,
    Column 6: "strand"
      The strand of the feature. + for positive strand , - for minus strand, and . for features that
      are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown.
    Column 7: "phase"
      For features of type "CDS", the phase indicates where the next codon begins relative to the 5'
      end of the current CDS feature. The 5' end for CDS features on the plus strand is the feature
      start and the 5' end for CDS features on the minus strand is the feature's end. The phase is
      0, 1, or 2, indicating the number of bases forward from the start of the current CDS feature
      the next codon begins. A phase of "0" indicates that a codon begins on the first nucleotide of
      the CDS feature (i.e. 0 bases forward). Note that ‘Phase’ in the context of a GFF3 CDS feature
      should not be confused with the similar concept of frame that is also a common concept in
      bioinformatics. Frame is generally calculated as a value for a given base relative to the
      start of the complete open reading frame (ORF).

      The phase is REQUIRED for all CDS features.
    Column 8: GTF/GFF2
     A semicolon-separated list of tag-value pairs, providing additional information about each
     feature. Except for the terminal semicolon, semicolons should be followed by spaces

     Each attribute is a key:value pair. Keys are expected to be single tokens; Values must be
     enclosed in parentheses so that spaces may be present. Attributes may contain leading/trailing
     spaces. gene_id and transcript_id are required, but may be empty strings.

    Column 8: GFF3
      A semicolon-separated list of feature attributes in the format tag=value. Multiple tag=value
      pairs are separated by semicolons. URL escaping rules are used for tags or values containing
      the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced
      with the %09 URL escape. Attribute values do not need to be and should not be quoted unless
      they are actually part of the value. Quotes should be included as part of the value by parsers
      and not stripped.

      For attributes with predefined meanings see
      https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

    ============================================================================================="""
    # Predefined keys for columns 0-8, agrees with GFF3
    column = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute']

    def __init__(self, row, format="gtf"):

        self.parsed = ''
        self.del_attr = None

        # attr_sep separates the key/value pairs in attribute
        self.attr_sep = '='
        if format == 'gtf' or 'gff2':
            self.attr_sep = ' '

    def feature_parse(self):
        """-----------------------------------------------------------------------------------------
        parse a feature line. Attributes are split on semicolons, and key/value pairs split on
        attr_sep character.

        :return:
        -----------------------------------------------------------------------------------------"""
        field = self.row.split(maxsplit=8)
        parsed = {}
        for i in 9:
            # extract the 9 defined columns
            if i in [3, 4]:
                # change begin, end, and score to int
                # should score be float?
                parsed[GxfRecord.column[i]] = int(field[i])
            else:
                parsed[GxfRecord.column[i]] = field[i]

        # split the attributes on ; and save as a hash

        field = parsed['attribute'].rstrip().split(';')
        # attribute may end in; so last field may be blank
        if not field[-1]:
            field.pop()

        for f in field:
            (key, value) = f.strip().split(self.attr_sep, maxsplit=1)
            parsed[key] = value.strip('"')

        if self.del_attr:
            del parsed['attribute']

        self.parsed = parsed
        return parsed

    def comment_parse(self):
        """-----------------------------------------------------------------------------------------
        intended for parsing comments; not implemented
        :return:
        -----------------------------------------------------------------------------------------"""
        pass
        return True


class GxfSet:
    """=============================================================================================
    A container for GxfRecords with tools for selecting sets of features from files or GxfSets

    ============================================================================================="""


# ##################################################################################################
# Testing
# #################################################################################################
if __name__  ## '__main__':
    exit(0)

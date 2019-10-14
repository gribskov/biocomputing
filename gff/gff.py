class Gff:
    ####################################################################################################
    # gff.py
    #
    # 9 October 2019    Michael Gribskov
    ####################################################################################################
    import sys
    import re

    column = ['sequence', 'method', 'feature', 'begin', 'end', 'score', 'strand', 'frame',
              'attribute']

    def __init__(self, file=""):
        self.data = []
        self.gff_in = None
        attr_sep = ' '

        if file:
            fh = self.open(file)
            if fh:
                self.gff_in = fh

    def open(self, file):
        """-----------------------------------------------------------------------------------------
        open a file for reading

        :param self:
        :param file: str, path to a GFF file
        :return: fh or False if unsuccessful
        -----------------------------------------------------------------------------------------"""
        try:
            fh = open(file, 'r')
            return fh
        except IOError:
            sys.stderr.write("gff.open unable to open GFF file ({})".format(file))
            return False

    def read(self):
        """-----------------------------------------------------------------------------------------
        Read a line from gff_in and delegate to the proper parsing function

        :return:
        -----------------------------------------------------------------------------------------"""
        self.line = self.gff_in.readline().rstrip()

        if self.line:
            if self.line.startswith('#'):
                self.comment_parse()
            else:
                self.feature_parse()

            return True

        else:
            # EOF
            return False

    def feature_parse(self):
        """-----------------------------------------------------------------------------------------
        parse a feature line
        :return:
        -----------------------------------------------------------------------------------------"""
        field = self.line.split(maxsplit=8)
        parsed = {}
        for i in range(len(field)):
            # extract the 9 defined columns
            parsed[Gff.column[i]] = field[i]

        # split the attributs on ; and restore as a hash

        field = parsed['attribute'].rstrip().split(';')
        # attribute may end in; so last field may be blank
        if not field[-1]:
            field.pop()

        for f in field:
            (key, value) = f.strip().split(self.attr_sep, maxsplit=1)
            parsed[key] = value.replace('"', '')

        self.data.append(parsed)

        return

    def comment_parse(self):
        pass
        return True

    def get_by_feature(self, key, start=0):
        """-----------------------------------------------------------------------------------------
        Generator for lines that match a feature tag

        :param key: str, string matching a feature in column 3
        :param start: int, line on which to start
        :return:
        -----------------------------------------------------------------------------------------"""
        data = self.data
        for n in range(start, len(data)):
            if data[n]['feature'] == key:
                yield n, data[n]

        raise StopIteration

    def get_by_value(self, column, key, start=0, stop=0):
        """-----------------------------------------------------------------------------------------
        A generator that returns rows where the speicified column matches the specified value.
        Rows missing columns, e.g., those generated from attributes, are skipped

        :param column: str, predefined or attribute column
        :param key: str, column value to match
        :param start: int, beginning row
        :param stop: int, ending row + 1
        :return: row, content_hash
        -----------------------------------------------------------------------------------------"""
        data = self.data
        if stop == 0:
            stop = len(data)

        for n in range(start, stop):
            if column not in data[n]:
                continue

            if data[n][column] == key:
                yield n, data[n]

        raise StopIteration

    def replace_by_column(self, column, find, replace):
        """-----------------------------------------------------------------------------------------
        Replace all of find by replace in column

        :param column: str, name of a predifined or attribute column
        :param find: str, string to replace
        :param replace: string to substitute for find
        :return: int, rows examined
        -----------------------------------------------------------------------------------------"""
        data = self.data
        n = 0
        for d in data:
            if column in d:
                d[column] = d[column].replace(find, replace)
                n += 1

        return n

    def replace_columns_re(self, column_list, search, replace=''):
        """-----------------------------------------------------------------------------------------
        replace strings in a list of columns using a regular expression
        -----------------------------------------------------------------------------------------"""
        data = self.data
        n = 0

        query = re.compile(search)
        for d in data:
            for column in column_list:
                if column in d:
                    d[column] = query.sub(replace, d[column])
                    n += 1

        return n

# ==================================================================================================
# test
# ==================================================================================================
if __name__ == '__main__':
    import sys
    import re

    # gff = Gff(file='stringtie.gff')
    gff = Gff(file='genome.gff')
    gff.attr_sep = '='

    line = 0
    while gff.read():
        line += 1

    sys.stdout.write('{} lines read\n'.format(line))

    # remove the string 'lcl|' in the sequence names
    # gff.replace_by_column('sequence', 'lcl|', '')
    query = r'(maker|augustus|masked|processed|gene)(-|_)'
    gff.replace_columns_re(['Parent', 'ID'], query, r'')

    # transcripts is a generator function
    # transcripts = gff.get_by_feature('transcript')
    transcripts = gff.get_by_feature('mRNA')
    (begin, transcript) = next(transcripts)

    for (end, transcript) in transcripts:
        # sys.stdout.write('{}\t{}\t{}\n'.format(transcript['gene_id'],
        #                                        transcript['transcript_id'],
        #                                        transcript['sequence']
        #                                        ))
        sys.stdout.write('{}\t{}\t{}\n'.format(transcript['Parent'],
                                               transcript['ID'],
                                               transcript['sequence']
                                               ))

        for (exon_n, exon) in gff.get_by_value('feature', 'exon', begin, end):
            print('\t{}\t{}\t{}\t{}'.
                  # format(exon['exon_number'], exon['begin'], exon['end'], exon['strand']))
                  format(exon['ID'], exon['begin'], exon['end'], exon['strand']))

        begin = end

    # the final transcript
    sys.stdout.write('{}\t{}\t{}\n'.format(transcript['Parent'],
                                           transcript['ID'],
                                           transcript['sequence']
                                           ))
    for (exon_n, exon) in gff.get_by_value('feature', 'exon', begin, len(gff.data)):
        print('\t{}\t{}\t{}\t{}'.
              # format(exon['exon_number'], exon['begin'], exon['end'], exon['strand']))
              format(exon['ID'], exon['begin'], exon['end'], exon['strand']))

    # print(flist)

    exit(0)

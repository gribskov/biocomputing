class Gff:
    ####################################################################################################
    # gff.py
    #
    # 9 October 2019    Michael Gribskov
    ####################################################################################################
    import sys

    column = ['sequence', 'method', 'feature', 'begin', 'end', 'score', 'strand', 'frame',
              'attribute']

    def __init__(self, file=""):
        self.data = []
        self.gff_in = None

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
        # attribute ends in; so last field is blank
        field.pop()
        for f in field:
            (key, value) = f.strip().split(' ', maxsplit=1)
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


# ==================================================================================================
# test
# ==================================================================================================
if __name__ == '__main__':
    import sys

    gff = Gff(file='stringtie.gff')

    line = 0
    while gff.read():
        line += 1

    sys.stdout.write('{} lines read\n'.format(line))

    flist = list(gff.get_by_feature('transcript'))
    (tnum, transcript) = flist[0]
    begin = tnum

    for tnum in range(1, len(flist)):
        sys.stdout.write('{}\t{}\n'.format(transcript['gene_id'], transcript['transcript_id']))
        (end, transcript) = flist[tnum]

        for (exon_n,exon) in gff.get_by_value('feature', 'exon', begin, end):
            print('\t{}\t{}\t{}\t{}'.
                  format(exon['exon_number'], exon['begin'], exon['end'], exon['strand']))

        begin = end

    sys.stdout.write('{}\t{}\n'.format(transcript['gene_id'], transcript['transcript_id']))
    for (exon_n,exon) in gff.get_by_value('feature', 'exon', begin, len(gff.data)):
        print('\t{}\t{}\t{}\t{}'.
                  format(exon['exon_number'], exon['begin'], exon['end'], exon['strand']))

    # print(flist)

    exit(0)

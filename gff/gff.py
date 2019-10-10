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
        self.line = self.gff_in.readline()

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
        attr = {}
        for f in field:
            (key, value) = f.strip().split(' ', maxsplit=1)
            attr[key] = value.replace('"', '')

        parsed['attribute'] = attr

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


# ==================================================================================================
# test
# ==================================================================================================
if __name__ == '__main__':
    import sys

    gff = Gff(file='stringtie.gff')

    line = 0
    while gff.read():
        line += 1

    sys.stdout.write('{} lines read'.format(line))

    for (line, entry) in gff.get_by_feature('transcript'):
        print(line, entry)

    exit(0)

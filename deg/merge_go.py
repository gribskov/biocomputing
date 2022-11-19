"""=================================================================================================
gather gene sets and carry over to genes identified by blast matching. Inputs
1. One or more GO source with source ID and lists of GO terms
2. blast search mapping source IDs  mapping source ID to new ID
3. gene ontology in OBO format to extend leaf GO terms to complete hierarchy

Michael Gribskov     18 November 2022
================================================================================================="""
import sys


class GOset():
    """=============================================================================================

    ============================================================================================="""

    def __init__(self, gofile=''):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.source = ''
        self.description = ''
        self.fh = None
        self.nlines = 0
        self.term = {}
        self.tag = {'GO': 'GO:'}
        self.column = []

        if gofile:
            self.source = gofile
            self.openfile(gofile, 'r')

    def openfile(self, filename, mode='', description=''):
        """-----------------------------------------------------------------------------------------
        safely open file and set source, description and fh

        :param filename: string     file available for reading
        :param description: string  text description of file content
        :return: file handle
        -----------------------------------------------------------------------------------------"""
        try:
            fh = open(filename, mode)
        except OSError:
            sys.stderr.write(f'GOset::openfile cannot open file ({filename})\n')
            exit(1)

        self.fh = fh
        self.description = description
        self.nlines = 0

        return fh

    def read_info(self, line, tagsplit='|'):

        """---------------------------------------------------------------------------------------------
        MSU     Soltu.DM.02G022630.1    Soltu.DM.02G022630.1            GO:0000785      TAIR:AT1G20696  IEA             C                       gene    taxon:4113      20201112        MSU
        MSU     Soltu.DM.02G022630.1    Soltu.DM.02G022630.1            GO:0003682      TAIR:AT1G20696  IEA             F                       gene    taxon:4113      20201112        MSU
        MSU     Soltu.DM.02G022630.1    Soltu.DM.02G022630.1            GO:0030527      TAIR:AT1G20696  IEA             F                       gene    taxon:4113      20201112        MSU

        Soltu.DM.01G045700.1    0e68191f38931e266e404adc8d317738        552     TIGRFAM TIGR01035       hemA: glutamyl-tRNA reductase   101     509     4.0E-114        T       11-11-2020      IPR000343     Tetrapyrrole biosynthesis, glutamyl-tRNA reductase       GO:0008883|GO:0033014|GO:0050661|GO:0055114

        :return:
        -----------------------------------------------------------------------------------------"""
        info = {}
        column = self.column

        if line:
            self.nlines += 1
            field = line.rstrip().split()
            for col in column:
                info[col] = field[column[col]]

            for tag in self.tag:
                # check for all tags
                target = self.tag[tag]
                tagset = []

                for f in field:
                    # check all columns for tagged data
                    if f.find(target) > -1:
                        if f.find(tagsplit):
                            tagset += f.split(tagsplit)
                        else:
                            tagset += f

                info[tag] = tagset

        self.add_info(info)

        return

    def add_info(self, info):
        """-----------------------------------------------------------------------------------------
        merge information about the current GO term into the existing information

        :param info: dict       source_id, tair, GO:
        :return: int            number of GO terms for this source_id
        -----------------------------------------------------------------------------------------"""
        try:
            id = info['source_id']
        except KeyError:
            return 0

        if not id in self.term:
            # unknown source term, create columns
            self.term[id] = {x: None for x in self.column}
            for t in self.tag:
                self.term[id][t] = []

        entry = self.term[id]
        for col in info:
            if col == 'GO':
                entry[col] += info['GO']
            else:
                entry[col] = info[col]

        return len(entry['GO'])

    def load(self, comment='#'):
        """-----------------------------------------------------------------------------------------
        load entire data file
        :return:
        -----------------------------------------------------------------------------------------"""
        self.nlines = 0
        for line in self.fh:
            if line.startswith(comment):
                continue

            self.read_info(line)

        return self.nlines

    def golist(self):
        """-----------------------------------------------------------------------------------------
        Return a unique list of GO terms in stored for all source IDs

        :return: list   unique GO terms
        -----------------------------------------------------------------------------------------"""
        go = set()
        for source in self.term:
            go.update(self.term[source]['GO'])

        return sorted(go)

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gofile = sys.argv[1]
    go = GOset(gofile)
    go.column = {'source_id': 1, 'tair': 4}
    tags = ['GO:']
    nlines = go.load()
    sys.stderr.write(f'processed {nlines} lines from {gofile}\n')
    sys.stderr.write(f'{len(go.term)} sequence IDs\n')

    golist = go.golist()
    sys.stderr.write(f'{len(golist)} unique GO terms\n')

    exit(0)

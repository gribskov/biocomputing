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
    # read GO file based on blast comparison to TAIR
    gofile = sys.argv[1]
    sys.stderr.write(f'reading GO assignments from {gofile}\n')

    go = GOset(gofile)
    go.column = {'source_id': 1, 'tair': 4}
    nlines = go.load()
    sys.stderr.write(f'processed {nlines} lines from {gofile}\n')
    sys.stderr.write(f'{len(go.term)} sequence IDs\n')
    golist = go.golist()
    sys.stderr.write(f'{len(golist)} unique GO terms\n')

    # read go file based on interproscan
    gofile = sys.argv[2]
    sys.stderr.write(f'\nreading GO assignments from {gofile}\n')
    goipr = GOset(gofile)
    goipr.column = {'source_id': 0}
    nlines = goipr.load()
    sys.stderr.write(f'processed {nlines} lines from {gofile}\n')
    sys.stderr.write(f'{len(goipr.term)} sequence IDs\n')
    golist = goipr.golist()
    sys.stderr.write(f'{len(golist)} unique GO terms\n')

    # read in list of trinity names
    trinity = []
    seqfile = '../sequential/filtered.names.txt'
    seq = open(seqfile, 'r')
    for line in seq:
        if not line:
            continue
        trinity.append(line.rstrip())
    seq.close()

    sys.stderr.write(f'\n{len(trinity)} sequences read from {seqfile}\n')

    # read mapping between DM_1-3_516_R44_potato.v6.1 and trinity from blastx search
    blastfile = 'best_blast_DM_1-3_516_R44_potato.v6.1.working_models.pep_1e-5.dmndblastx'
    blast = open(blastfile, 'r')
    dm2trinity = {}
    blastn = 0
    trinityn = 0
    dupn = 0
    for line in blast:
        blastn += 1
        field = line.split()
        if field[0] not in trinity:
            continue

        trinityn += 1
        if field[4] in dm2trinity:
            dm2trinity[field[4]].append(field[0])
            dupn += 1
            sys.stderr.write(f'\t{field[4]} is duplicate {field[0]}\n')
        else:
            dm2trinity[field[4]] = [field[0]]

    sys.stderr.write(f'{blastn} blast results processed\n')
    sys.stderr.write(f'{trinityn} trinity genes found\n')
    sys.stderr.write(f'{dupn} duplicates found\n')
    sys.stderr.write(f'{len(dm2trinity)} unique dm genes\n')

    # construct a new gene to ontology mapping indexed by trinity name
    tgo = GOset()
    trinity2dm = {}
    for gene in trinity:
        tgo.term[gene] = {'GO': []}
        trinity2dm = []

    for dm in dm2trinity:
        for trinity in dm2trinity[dm]:
            if dm in go.term:
                golist = set(go.term[dm]['GO'])
            if dm in goipr.term:
                golist.update(goipr.term[dm]['GO'])
        tgo.term[trinity]['GO'] = list(golist)

    # make sure all terms include the biological process, molecular function and cellular component roots
    # so that the terms will never show up as significant
    for trinity in tgo.term:
        for go in ('GO:0008150', 'GO:0003674', 'GO:0005575'):
            if go not in tgo.term[trinity]['GO']:
                tgo.term[trinity]['GO'].append(go)

    # write out as a read mapping: tab delimited gene followed by list of GO terms
    # tab = '\t'
    # mergeout = open('merged_go', 'w')
    # for trinity in tgo.term:
    #     mergeout.write(f'{trinity}')
    #     for go in tgo.term[trinity]['GO']:
    #         mergeout.write(f'{tab}{go}')
    #     mergeout.write('\n')
    # mergeout.close()

    tab = '\t'
    mergeout = open('merged_go', 'w')
    for trinity in tgo.term:
        for go in tgo.term[trinity]['GO']:
            mergeout.write(f'{trinity}{tab}{go}\n')

    mergeout.close()

    exit(0)

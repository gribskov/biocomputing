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

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.source = ''
        self.description = ''
        self.fh = None
        self.nlines = 0
        self.term = {}

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

    def read_info(self, column={}, tag={}, tagsplit='|'):

        """---------------------------------------------------------------------------------------------
        MSU     Soltu.DM.02G022630.1    Soltu.DM.02G022630.1            GO:0000785      TAIR:AT1G20696  IEA             C                       gene    taxon:4113      20201112        MSU
        MSU     Soltu.DM.02G022630.1    Soltu.DM.02G022630.1            GO:0003682      TAIR:AT1G20696  IEA             F                       gene    taxon:4113      20201112        MSU
        MSU     Soltu.DM.02G022630.1    Soltu.DM.02G022630.1            GO:0030527      TAIR:AT1G20696  IEA             F                       gene    taxon:4113      20201112        MSU

        Soltu.DM.01G045700.1    0e68191f38931e266e404adc8d317738        552     TIGRFAM TIGR01035       hemA: glutamyl-tRNA reductase   101     509     4.0E-114        T       11-11-2020      IPR000343     Tetrapyrrole biosynthesis, glutamyl-tRNA reductase       GO:0008883|GO:0033014|GO:0050661|GO:0055114

        :return:
        ---------------------------------------------------------------------------------------------"""
        line = self.fh.readline()
        info = {}

        if line:
            self.nlines += 1
            field = line.rstrip.split()
            for col in column:
                info[col] = field[column[col]]

            for tag in tags:
                # check for all tags
                tagset = []

                for f in field:
                    # check all columns for tagged data
                    if f.find(tag) > -1:
                        if f.find(tagsplit):
                            tagset += f.split(tagsplit)
                        else:
                            tagset += f

                info[tag] = tagset

        self.addinfo(info)


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    exit(0)

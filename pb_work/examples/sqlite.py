"""=================================================================================================
SQLite3 example.  This example combines a GFF object and the same database as in sqlite2.py

Michael Gribskov     24 February 2019
================================================================================================="""
import sys
import sqlite3 as sq3


class Genome():
    """---------------------------------------------------------------------------------------------
    Persistent genome class, based on SQLite3
    ---------------------------------------------------------------------------------------------"""

    def __init__(self, dbfile='genome.db', new=False):
        """-----------------------------------------------------------------------------------------
        creates the database connection and initializes tables

        Database is set up to use the SQLite row interface
        -----------------------------------------------------------------------------------------"""
        self.dbh = sq3.connect(dbfile)
        self.db = self.dbh.cursor()
        self.db.row_factory = sq3.Row

        self.setupDBTables(new)

        # end of __init__

    def setupDBTables(self, new):
        """-----------------------------------------------------------------------------------------
        Creates all tables for the database
        :param db: database cursos from SQLite3
        :return: None
        -----------------------------------------------------------------------------------------"""
        if new:
            self.db.execute('DROP TABLE IF EXISTS chromosome')
            self.db.execute('DROP TABLE IF EXISTS gene')
            self.db.execute('DROP TABLE IF EXISTS exon')

        sql = '''
             CREATE TABLE IF NOT EXISTS chromosome
             (   chr_id TEXT PRIMARY KEY,
                 length INTEGER
             )   
    
             '''
        self.db.execute(sql)

        sql = '''
             CREATE TABLE IF NOT EXISTS gene
             (   gene_id TEXT PRIMARY KEY,
                 chr_id TEXT,
                 begin INTEGER,
                 end INTEGER,
                 strand TEXT
             )   
    
             '''
        self.db.execute(sql)

        sql = '''
             CREATE TABLE IF NOT EXISTS exon
             (   exon_id TEXT PRIMARY KEY,
                 gene_id TEXT,
                 begin INTEGER,
                 end INTEGER,
                 strand TEXT,
                 type
             )   
    
             '''
        self.db.execute(sql)

        return None


class Gff:
    """=============================================================================================
    Simple class to read features from GFF3 file and store in database

    Michael Gribskov    24 February 2019
    ============================================================================================="""

    def __init__(self, filename, dbfilename='genomedb.db', dbnew=True):
        """-----------------------------------------------------------------------------------------
        Gff constructor

        :param filename: string, GFF3 file name
        -----------------------------------------------------------------------------------------"""
        self.filename = filename
        self.fh = None

        self.dbfilename = dbfilename
        self.dbnew = dbnew
        self.genome = Genome(self.dbfilename, dbnew)

        try:
            self.fh = open(self.filename, 'r')
            self.Gff_load()

        except IOError:
            sys.stderr.write('Gff::__init__ - error opening file ({})\n'.format(self.filename))

    def Gff_load(self):
        """-----------------------------------------------------------------------------------------
        load a GFF3 data file into the database
        GFF3 format is defined at https://useast.ensembl.org/info/website/upload/gff3.html

        :return: int, number of rows loaded
        -----------------------------------------------------------------------------------------"""
        nrow = 0
        for line in self.fh:
            if line.startswith('#'):
                continue

            nrow += 1
            (seqid, source, type, start, end, score, strand, phase, attributes) = \
                line.rstrip().split('\t')
            # except ValueError:
            #     print('Error', line)
            #     continue

            if type.find('gene') >= 0:
                # load gene into gene table
                attr = self.parse_attributes(attributes)
                sql = '''
                INSERT INTO gene
                VALUES ( "{}", "{}", {}, {}, "{}" );
                '''.format(attr['ID'], seqid, start, end, strand)

                self.genome.db.execute(sql)

        return nrow

    def parse_attributes(self, attributes):
        """-----------------------------------------------------------------------------------------
        Attributes are a semicolon delimited string, with each clause a tag/value pair separated
        by equals

        :param attributes: string
        :return: dict, parsed attributes
        -----------------------------------------------------------------------------------------"""
        field = attributes.split(';')
        attr = {}
        for f in field:
            tag, value = f.split('=')
            attr[tag] = value

        return attr


# --------------------------------------------------------------------------------------------------
# main/testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gff = Gff('../2019/HW5/at_1000k.gff3')

    print(gff)

    exit(0)

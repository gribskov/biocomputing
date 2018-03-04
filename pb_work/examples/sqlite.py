"""=================================================================================================
SQLite3 example

Michael Gribskov     March 04  2018
================================================================================================="""
import sqlite3 as sq3


class Genome():
    """---------------------------------------------------------------------------------------------
    Persistent genome class, based on SQLite3
    ---------------------------------------------------------------------------------------------"""

    def __init__(self, dbfile='genome.db', new=False):
        """-----------------------------------------------------------------------------------------
        creates the database connection and initializes tables
        -----------------------------------------------------------------------------------------"""
        self.dbh = sq3.connect(dbfile)
        self.db = self.dbh.cursor()

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


# --------------------------------------------------------------------------------------------------
# main/testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    genome = Genome('mydbfile.db', new=True)

    exit(0)

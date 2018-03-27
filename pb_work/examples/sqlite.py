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

    chromosomes = [ ('1', 30427671),
                    ('2', 19698289),
                    ('3', 23459830),
                    ('4', 18585056),
                    ('5', 26975502)
                   ]

    genome.db.executemany('INSERT INTO chromosome VALUES (?, ?)', chromosomes)
    genome.dbh.commit()

    sql = '''
    SELECT *
    FROM chromosome
    '''
    genome.db.execute(sql)
    for row in genome.db.fetchall():
        print(row[0], row[1])

    genome.db.row_factory = sq3.Row
    genome.db.execute(sql)
    for row in genome.db:
        for key in row.keys():
            print('{}:{}'.format(key,row[key]),end='\t')
        print()

    exit(0)

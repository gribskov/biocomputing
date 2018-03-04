"""=================================================================================================
SQLite3 example

Michael Gribskov     March 04  2018
================================================================================================="""
import sqlite3 as sq3


def setupDBTables(db):
    """---------------------------------------------------------------------------------------------
    Creates all tables for the database
    :param db: database cursos from SQLite3
    :return: None
    ---------------------------------------------------------------------------------------------"""
    sql = '''
         CREATE TABLE IF NOT EXISTS chromosome
         (   chr_id TEXT PRIMARY KEY,
             length INTEGER
         )   

         '''
    db.execute(sql)

    sql = '''
         CREATE TABLE IF NOT EXISTS gene
         (   gene_id TEXT PRIMARY KEY,
             chr_id TEXT,
             begin INTEGER,
             end INTEGER,
             strand TEXT
         )   

         '''
    db.execute(sql)

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
    db.execute(sql)

    sql = '''
         INSERT INTO chromosome 
             VALUES (1, 30427671 )
         '''
    db.execute(sql)

    return None


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    dbh = sq3.connect('mydbfile.db')
    db = dbh.cursor()
    setupDBTables(db)

    dbh.commit()

    exit(0)

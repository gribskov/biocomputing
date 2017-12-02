'''--------------------------------------------------------------------------------
kollemadb
Database object for kollema transcriptome annotation tool
--------------------------------------------------------------------------------'''
import sys
import sqlite3 as sql


class Kollemadb(object):
    def __init__(self, dbfile=':memory:', new=False):
        '''

        :param dbfile: database file, use :memory: for in memory only
        :param new: initialize tables (be careful)
        :return: always true for now
        '''
        # isolation_level None is autocommit
        dbh = sql.connect(dbfile, isolation_level=None)
        dbh.row_factory = sql.Row
        self.db = dbh.cursor()

        if new:
            self.initAllTables()

        return None

    def initAllTables(self):
        '''

        :return:
        '''
        self.initProjectTable()

        return True

    def initProjectTable(self):
        '''

        :return:
        '''
        sql = '''
            DROP TABLE IF EXISTS project;
            CREATE TABLE project (
                name TEXT,
                owner TEXT,
                created TEXT,
                status TEXT,
                updated TEXT );
        '''
        self.db.executescript(sql)

        return True

'''--------------------------------------------------------------------------------
kollemadb testing

--------------------------------------------------------------------------------'''
if __name__ == '__main__':

    kdb = Kollemadb(new=True)

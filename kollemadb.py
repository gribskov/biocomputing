'''--------------------------------------------------------------------------------
kollemadb
Database object for kollema transcriptome annotation tool
--------------------------------------------------------------------------------'''
import sys
import re
import sqlite3 as s3


class Kollemadb(object):
    def __init__(self, dbfile=':memory:', new=False):
        '''

        :param dbfile: database file, use :memory: for in memory only
        :param new: initialize tables (be careful)
        :return: always true for now
        '''
        # isolation_level None is autocommit
        dbh = s3.connect(dbfile, isolation_level=None)
        dbh.row_factory = s3.Row
        self.db = dbh.cursor()

        if new:
            self.initAllTables()

        return None

    def initAllTables(self):
        '''
        master build called to initialize all tables
        :return:
        '''
        self.initProjectTable()
        self.initTaskTable()

        return True

    def initProjectTable(self):
        '''
        sql to build the project table
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

    def initTaskTable(self):
        '''
        sql to build the task table
        '''
        sql = '''
            DROP TABLE IF EXISTS task;
            CREATE TABLE task (
                projectid INTEGER,
                name TEXT,
                begin TEXT,
                end TEXT,
                filename TEXT
                );
        '''
        self.db.executescript(sql)

        return True

    def showTables(self):
        '''
        list the tables in the database and their column definition
        usage
            print( kollemadb.showTables() )
        '''
        rtab = re.compile('\t$')
        rspace = re.compile('\s+')
        s = ''
        self.db.execute("SELECT * FROM sqlite_master")
        #self.db.execute("SELECT * FROM sqlite_master WHERE type='table'")
        for row in kdb.db:
            for r in row:
                if isinstance(r, str):
                    r = r.strip()
                    r = rspace.sub(' ',r)
                s += '{0}\t'.format(r)
            #s+= '\n'
            s = rtab.sub('\n',s)
        return s.strip()

    def get(self, table, limits=''):
        '''
        retrieve contents of table as a dictionary
        usage
            project = kollemadb.get(project)
            task = kollemadb.get(task)
        '''
        sql = '''
            SELECT * FROM {0}
            '''.format(table)
        
        try:
            self.db.execute(sql)
        except s3.Error as e:
            print("Kollemadb::get - error:", e.args[0])
            #sqlite3.OperationalError: no such table:

        result = {}
        row = self.db.fetchone()

        for k in row.keys():
            result[k] = row[k]

        return result


    def set(self, table, data):
        '''
        add a new row to an existing table
        usage
            kollemadb.set( table, data )
        '''
        columns = ''
        values = ''
        comma = ''
        for k in data:
            columns += '{0}{1} '.format(comma, k)
            values += "{0}'{1}' ".format(comma, data[k])
            comma = ','

        sql = '''
            INSERT INTO {0} 
                ({1}) VALUES ({2})'''.format(table, columns.strip(), values.strip())

        print('sql:', sql)
        try:
            self.db.execute(sql)
        except s3.Error as e:
            print("Kollemadb::set - error:", e.args[0])



'''--------------------------------------------------------------------------------
kollemadb testing

--------------------------------------------------------------------------------'''
if __name__ == '__main__':

    kdb = Kollemadb(new=True)
    print(kdb.showTables())

    data = { 'name':'test',
             'owner':'gribskov',
             'status':'old'
           }
    kdb.set('project', data)
    result = kdb.get('project')
    for k in result:
        print(k, result[k])

    #result = kdb.get('task')
    #result = b.get('transcript')

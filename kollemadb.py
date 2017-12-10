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

        :dbfile: database file, use :memory: for in memory only, default = :memory:
        :new: initialize tables (be careful), default = False
        :return: None
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
        :return: True
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
        # self.db.execute("SELECT * FROM sqlite_master WHERE type='table'")
        for row in kdb.db:
            for r in row:
                if isinstance(r, str):
                    r = r.strip()
                    r = rspace.sub(' ', r)
                s += '{0}\t'.format(r)
            # s+= '\n'
            s = rtab.sub('\n', s)
        return s.strip()

    def get(self, table, limits=''):
        '''
        retrieve contents of table as a dictionary
        :param table: table name
        :param limits: sql where clause
        :return: number of rows returned by query
        usage
            project = kollemadb.get('project')
            task = kollemadb.get('task')
        '''
        sql = '''
            SELECT * FROM {0}
            '''.format(table)

        try:
            self.db.execute(sql)
        except s3.Error as e:
            print("Kollemadb::get - error:", e.args[0])
            # sqlite3.OperationalError: no such table:

        result = []
        for row in self.db.fetchall():
            r = {}
            for k in row.keys():
                r[k] = row[k]
            result.append(r)

        # store result in object
        self.result = result

        return len(result)

    def asFormatted(self, indent=2):
        '''

        :param table: table name
        :param limits: SQL limits (where clause)
        :return:
        '''
        out = ''
        if len(self.result) == 0:
            return out

        # get the length of the longest entry in each column
        pad = 2
        fieldsize = {}
        for row in self.result:
            totallen = 0
            for k in row.keys():
                if not k in fieldsize:
                    fieldsize[k] = len(k)
                if row[k] == None:
                    continue
                if fieldsize[k] < len(row[k]):
                    fieldsize[k] = len(row[k])
                totallen += fieldsize[k] + pad + 1
        totallen += 1

        nrow = 0
        # print('{}'.format('-'*totallen))
        out += '{}{}\n{}'.format(' ' * indent, '-' * totallen, ' ' * indent)
        for row in self.result:
            nrow += 1
            out += '|'
            if nrow == 1:
                # print column titles
                for k in row.keys():
                    f = '{:<' + '{}'.format(fieldsize[k] + pad) + '}|'
                    # print('f=', f)
                    out += f.format(k)

                # print(str)
                out += '\n'
                out += '{}{}\n{}'.format(' ' * indent, '-' * totallen, ' ' * indent)
                out += '|'

            for k in row.keys():
                # print column values
                f = '{:<' + '{}'.format(fieldsize[k] + pad) + '}|'
                # print('f=', f)
                if row[k] == None:
                    out += f.format('')
                else:
                    out += f.format(row[k])
            out += '\n{}'.format(' ' * indent)
        out += '{}\n'.format( '-' * totallen )

        return out

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

        # print('sql:', sql)
        try:
            self.db.execute(sql)
        except s3.Error as e:
            print("Kollemadb::set - error:", e.args[0])

    def fromTerm(self, id):
        '''
        interactively acquire information for a table row
        :param id: table name
        :return: True
        '''
        self.db.execute('SELECT * FROM {0} LIMIT 1'.format(id))
        newdata = {}
        for tup in self.db.description:
            newdata[tup[0]] = input('{0} {1}:'.format(id, tup[0]))
        self.set(id, newdata)


'''--------------------------------------------------------------------------------
kollemadb testing

--------------------------------------------------------------------------------'''
if __name__ == '__main__':
    kdb = Kollemadb(new=True)
    print(kdb.showTables())

    # result = kdb.get('project')

    data = {'name': 'test',
            'owner': 'gribskov',
            'status': 'old'
            }
    kdb.set('project', data)

    kdb.fromTerm('project')
    kdb.get('project')
    kdb.asFormatted()

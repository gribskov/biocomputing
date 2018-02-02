"""---------------------------------------------------------------------------------------------------------------------
kollemadb
Database object for kollema transcriptome annotation tool
---------------------------------------------------------------------------------------------------------------------"""
import sys
import re
import sqlite3 as s3


class Kollemadb(object):
    """-----------------------------------------------------------------------------------------------------------------
    kollemadb
    Persistant class built on sqlite3
    -----------------------------------------------------------------------------------------------------------------"""

    def __init__(self, dbfile=':memory:', new=False):
        """
        Kollemadb class constructor
        :dbfile: database file, use :memory: for in memory only, default = :memory:
        :new: initialize tables (be careful), default = False
        :return: None
        """
        # isolation_level None is autocommit
        dbh = s3.connect(dbfile, isolation_level=None)
        dbh.row_factory = s3.Row
        self.db = dbh.cursor()

        if new:
            self.initAllTables()

    def initAllTables(self):
        """-------------------------------------------------------------------------------------------------------------
        Constructs all tables, does not drop existing tables if present
        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        self.initProjectTable()
        self.initTaskTable()
        self.initTranscriptTable()
        self.initPathTable()

        return True

    def initProjectTable(self):
        """-------------------------------------------------------------------------------------------------------------
        sql to build the project table
        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        sql = '''
            CREATE TABLE 
            IF NOT EXISTS
            project (
                name TEXT,
                owner TEXT,
                created TEXT,
                status TEXT,
                updated TEXT 
            );
        '''
        self.db.executescript(sql)

        return True

    def initTaskTable(self):
        """-------------------------------------------------------------------------------------------------------------
        sql to build the task table
        -------------------------------------------------------------------------------------------------------------"""
        sql = '''
            CREATE TABLE 
            IF NOT EXISTS
            task (
                projectid INTEGER,
                name TEXT,
                begin TEXT,
                end TEXT,
                filename TEXT
            );
        '''
        self.db.executescript(sql)

        return True

    def initTranscriptTable(self):
        """-------------------------------------------------------------------------------------------------------------
        sql to build transcript table
        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        sql = '''
            CREATE TABLE 
            IF NOT EXISTS
            transcript (
                id TEXT, 
                component INTEGER, 
                gene INTEGER, 
                isoform INTEGER,
                seq TEXT, 
                doc TEXT
            )
            '''
        self.db.execute(sql)

        return True

    def initPathTable(self):
        """-------------------------------------------------------------------------------------------------------------
        sql to build path table
        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        sql = '''
            CREATE TABLE 
            IF NOT EXISTS
            path (
                transcriptid INTEGER,
                step TEXT
            )
        
            '''
        self.db.execute(sql)

        return True

    def clear(self, table):
        """-------------------------------------------------------------------------------------------------------------
        delete all rows of the named table
        :param table: table name
        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        sql = '''
            DELETE FROM {table_id}
            '''.format(table_id=table)
        self.db.execute(sql)

        return True

    def clearAll(self):
        """-------------------------------------------------------------------------------------------------------------
        delete all rows from all tables
        :return: True
        -------------------------------------------------------------------------------------------------------------"""

        self.db.execute("SELECT * FROM sqlite_master WHERE type='table'")
        for row in self.db:
            self.clear(row[1])

        return True

    def loadFromList(self, table, columns, data):
        """
        Construct an SQL insert statement for the listed columns using placeholders
        use synchronize=False pragma and executemany to accelerate
        :param table: name of table
        :param columns: list of column names
        :param data: list of rows of values.  Each row is a list of data
        :return: number of rows
        """
        colnames = '({})'.format(','.join(columns))

        sql = '''
                INSERT INTO {id}
                {c}
                VALUES ({v})
                    '''.format(id=table, c=colnames, v=','.join('?' * len(columns)))
        res = self.db.executemany(sql, data)

        return res

    def showTables(self):
        """-------------------------------------------------------------------------------------------------------------
        list the tables in the database and their column definition
        usage
            print( kollemadb.showTables() )
        -------------------------------------------------------------------------------------------------------------"""
        rtab = re.compile('\t$')
        rspace = re.compile('\s+')
        s = ''
        self.db.execute("SELECT * FROM sqlite_master")
        # self.db.execute("SELECT * FROM sqlite_master WHERE type='table'")
        for row in self.db:
            for r in row:
                if isinstance(r, str):
                    r = r.strip()
                    r = rspace.sub(' ', r)
                s += '{0}\t'.format(r)
            # s+= '\n'
            s = rtab.sub('\n', s)
        return s.strip()

    def get(self, table, limits='', showsql=False):
        """-------------------------------------------------------------------------------------------------------------
        retrieve contents of table as a dictionary
        :param table: table name
        :param limits: sql where clause
        :param showsql: print out generated SQL if true
        :return: number of rows returned by query
        usage
            project = kollemadb.get('project')
            task = kollemadb.get('task')
        -------------------------------------------------------------------------------------------------------------"""
        sql = '''
                SELECT * FROM {0}
                {1}
                '''.format(table, limits)

        if showsql:
            print('kollemadb.get:{}'.format(sql))

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
        """-------------------------------------------------------------------------------------------------------------
        return string with a formatted tabular version of the current result
        usage
            print(kollemadb.asFormatted())
        :param indent: number of spaces to indent at the left
        :return: formatted string
        -------------------------------------------------------------------------------------------------------------"""
        out = ''
        if len(self.result) == 0:
            return out

        # get the length of the longest entry in each column
        pad = 2
        fieldsize = {}
        totallen = 0
        for row in self.result:
            totallen = 0
            for k in row.keys():
                if k not in fieldsize:
                    # initialize with column name
                    fieldsize[k] = len(k)
                if row[k] is None:
                    pass
                elif fieldsize[k] < len(str(row[k])):
                    fieldsize[k] = len(str(row[k]))
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
        out += '{}\n'.format('-' * totallen)

        return out

    def set(self, table, data):
        """-------------------------------------------------------------------------------------------------------------
        add a new row to an existing table
        usage
            kollemadb.set( table, data )
        -------------------------------------------------------------------------------------------------------------"""
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
        """-------------------------------------------------------------------------------------------------------------
        interactively acquire information for a table row
        :param id: table name
        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        self.db.execute('SELECT * FROM {0} LIMIT 1'.format(id))
        newdata = {}
        for tup in self.db.description:
            newdata[tup[0]] = input('{0} {1}:'.format(id, tup[0]))
        self.set(id, newdata)


'''---------------------------------------------------------------------------------------------------------------------
kollemadb testing

---------------------------------------------------------------------------------------------------------------------'''
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

    # test clearinng tables
    kdb.get('project')
    print(kdb.asFormatted())
    kdb.clearAll()
    kdb.get('project')
    print(kdb.asFormatted())

    cols = ['id', 'component', 'gene']
    data = [
        ('a', 1, 0),
        ('b', 1, 1),
        ('c', 1, 0)
    ]

    kdb.loadFromList('transcript', cols, data)
    kdb.get('transcript')
    print(kdb.asFormatted())

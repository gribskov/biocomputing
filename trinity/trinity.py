import re
from sequence.fasta import Fasta

lenre = re.compile('len=(\d+)')
pathre = re.compile('.*path=\[([^\]]*)\]')
idre = re.compile(r'>*TRINITY_([^_]+)_c(\d+)_g(\d+)_i(\d+)')


class Trinity(Fasta):
    """
    trinity transcript class.  Uses sqlite for persistance
    """

    # ----------------------------------------------------------------------------------------------------
    #  non-object utility functions
    # ----------------------------------------------------------------------------------------------------


    def getLen(line):
        """
        get the sequence length from the documentation. Documentation looks like
        len=230 path=[208:0-229] [-1, 208, -2]
        usage
            len = getLen(trinity.doc)
        """
        return lenre.match(line).group(1)

    def getPath(line):
        """
        return a list of path components from documentation. Documentation looks like
        len=230 path=[208:0-229] [-1, 208, -2]
        returns strings such as 208:0-229
        usage
            pathlist = getPath(trinity.doc)
        """
        path = pathre.match(line).group(1)
        plist = path.split(' ')
        return plist

    def splitID(line):
        """
        Breakdown the trinity ID string to give the
        Cluster,  component, gene and isoform
        usage
            cluster, component, gene, isoform = trinityID(trinity.id)
         """
        cluster, component, gene, isoform = idre.match(line).groups()
        return cluster, component, gene, isoform

    # ----------------------------------------------------------------------------------------------------
    #  # class methods
    # ----------------------------------------------------------------------------------------------------

    def __init__(self, dbfile=':memory:', new=True):
        super().__init__()
        self.dbname = ''
        self.db = ''
        self.rows = []
        self.clear()
        self.dbSetup(dbfile, new)

    def clear(self):
        '''
        clear the current transcript information. Does not affect sequence information
        usage
            self.clear()
        '''
        self.cluster = ''
        self.component = ''
        self.gene = ''
        self.isoform = ''
        self.path = []

    def dbSetup(self, dbfile=':memory:', new=True):
        """
        set up dataase connection.  Create tables if necessary.
        dbfile is either the database file name or ':memory:'
        usage
            trinity.dbSetup()
            trinity.dbSetup('trinity.sqlite')
            trinity.dbSetup('trinity.sqlite', new=False)
        """
        import sqlite3 as sql

        # isolation_level None is autocommit
        dbh = sql.connect(dbfile,isolation_level=None)
        dbh.row_factory = sql.Row
        self.db = dbh.cursor()

        if new:
            sql = '''
                DROP TABLE IF EXISTS transcript; 
                DROP TABLE IF EXISTS path
                '''
            self.db.executescript(sql)

            sql = '''
                CREATE TABLE 
                    IF NOT EXISTS
                    transcript (
                        cluster TEXT, 
                        component INTEGER, 
                        gene INTEGER, 
                        isoform INTEGER,
                        seq TEXT, 
                        doc TEXT
                    )
                '''
            self.db.execute(sql)

            sql = '''
                CREATE TABLE path 
                    (
                        transcriptid INTEGER,
                        step TEXT
                    )

                '''
            self.db.execute(sql)

    def insertTranscript(self):
        '''
        Insert the current transcript record into the database
        Affects transcript and path tables
        Returns rowid of new row in transcript table
        usage
            id = trinity.insertTranscript()
        '''
        sql_transcript = '''
        INSERT INTO transcript
        VALUES (?, ?, ?, ?, ?, ? )
        '''
        sql_path = '''
        INSERT INTO path
        VALUES (?,?)
        '''
        self.db.execute(sql_transcript, (self.cluster, self.component, self.gene, self.isoform, self.seq, self.doc))
        transcript_id = self.db.lastrowid

        for segment in self.path:
            self.db.execute(sql_path, (transcript_id, segment))

        return transcript_id

    def dumpTable(self, table):
        '''
        print out  the complete contents of the named table
        usage
        trinity.dumpTable('path')
        '''
        id_old = 0
        for row in self.db.execute('SELECT rowid, * FROM {0}'.format(table)):
            if not id_old == row[1]:
                print('{0} {1}'.format(table, row['rowid']))
                id_old = row[1]

            for key in row.keys():
                if key == 'rowid':
                    continue
                print('  {0}:{1}'.format(key, row[key]), end='')
            print()

    def formatFasta(self):
        '''
        return a string with the Fasta formatted version of the current transcript
        usage
            print(trinity.formatFasta)
        '''
        self.id = '{0}_c{1}_g{2}_i{3}'.format(self.cluster, self.component, self.gene, self.isoform)
        return self.format()

    def fromFasta(self, fastafile, limit=1000000000):
        '''
        populate the trinity datastructure from a file in Fasta format
        usage
        trinity.fromFasta( 'Trinity.fasta' )
        '''
        if not self.db:
            self.dbSetup()

        trinity.open(fastafile)
        ntranscript = 0;
        while trinity.next():
            trinity.clear()
            self.cluster, self.component, self.gene, self.isoform = Trinity.splitID(trinity.id)
            self.len = Trinity.getLen(trinity.doc)
            for segment in Trinity.getPath(trinity.doc):
                self.path.append(segment)
            row = self.insertTranscript()
            ntranscript += 1
            if ntranscript >= limit:
                break

        return row


if __name__ == '__main__':
    trinity = Trinity('trinity.sql')
    row = trinity.fromFasta(r'C:\Users\gribs\Dropbox\rice\Trinity.fasta', limit=205)
    print('{0} rows read from {1}'.format(row, trinity.filename))
    trinity.db.close()

    t2 = Trinity(dbfile='trinity.sql', new=False)
    t2.dumpTable('transcript')
    t2.dumpTable('path')

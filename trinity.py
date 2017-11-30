'''

'''
from sequence.fasta import Fasta
import re
import sqlite3 as sql

lenre = re.compile('len=(\d+)')
pathre = re.compile('.*path=\[([^\]]*)\]')


def getLen(line):
    '''
    get the sequence length from the documentation
    '''
    return lenre.match(line).group(1)


def getPath(line):
    '''
    return a list of path components from documentation
    '''
    # path = pathre.match(line).group(1)
    # print('doc:', line)
    # print('pre:',pathre)
    # print('path:', pathre.match(line).group(1))
    path = pathre.match(line).group(1)
    plist = path.split(' ')
    return plist


def trinityID(id):
    '''
    Breakdown the trinity ID string to give the
    Cluster,  component, gene and isoform
     '''
    cluster, component, gene, isoform = re.compile(r'>*TRINITY_([^_]+)_c(\d+)_g(\d+)_i(\d+)').match(id).groups()
    return cluster, component, gene, isoform


trinity = Fasta()

file = r'C:\Users\gribs\Dropbox\rice\Trinity.fasta'
trinity.fh = open(file, 'r')

# dbh = sql.connect('trinity.db')
dbh = sql.connect(":memory:")
dbh.row_factory = sql.Row
db = dbh.cursor()

sql = '''
DROP TABLE IF EXISTS transcript; DROP TABLE IF EXISTS path
'''
db.executescript(sql)

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
db.execute(sql)

sql = '''
CREATE TABLE path 
    (
        transcriptid INTEGER,
        step TEXT
    )

'''
db.execute(sql)

nseq = 0

while trinity.next():
    nseq += 1

    cluster, component, gene, isoform = trinityID(trinity.id)
    nn = getPath(trinity.doc)

    # print('cluster:', cluster, 'component:', component, 'gene:', gene, 'isoform:', isoform)
    db.execute('INSERT INTO transcript VALUES (?,?,?,?,?,?)',
               (cluster, component, gene, isoform, trinity.seq, trinity.doc))
    transcriptid = db.lastrowid
    for p in nn:
        db.execute('INSERT INTO path VALUES (?,?)', (transcriptid, p))
    if nseq > 10: break
    break

db.execute('SELECT id, component, gene, isoform FROM transcript LIMIT 5')
n = 0
for row in db:
    n += 1
    print(row['id'])
    for k in row.keys():
        print('    ', k, row[k])
        # print( '    ', row['cluster'])

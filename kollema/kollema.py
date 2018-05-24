"""---------------------------------------------------------------------------------------------------------------------
kollema
Transcriptome annotation
---------------------------------------------------------------------------------------------------------------------"""
import os
import sqlite3 as sql
from kollemadb import Kollemadb
from menu import Menu
from trinity.trinity import Trinity


def status():
    """-----------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------------------"""
    print('status activated')
    pass


def projectSelect(kdb, indent=2):
    """-----------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------------------"""
    project = input('\n{0}{1}: '.format(' ' * indent, 'Select project name'))
    response = kdb.get('project', limits='WHERE name = {}'.format(project))
    if response == 1:
        return project
    elif response == 0:
        print('Unknown project ({})'.format(project))
    else:
        print('Multiple projects called "{}" exist'.format(project))

    return ''


def transcriptLoad(kdb, indent=2):
    """-----------------------------------------------------------------------------------------------------------------
                id TEXT,
                component INTEGER,
                gene INTEGER,
                isoform INTEGER,
                seq TEXT,
                doc TEXT

    -----------------------------------------------------------------------------------------------------------------"""
    transcript_file = input('\n{0}{1}: '.format(' ' * indent, 'Transcript file to load'))
    trinity = Trinity()
    trinity.open(transcript_file)
    ntranscript = 0
    while trinity.next():
        kdb.set('transcript', {'id': trinity.cluster,
                               'component': trinity.component,
                               'gene': trinity.gene,
                               'isoform': trinity.isoform,
                               'seq': trinity.seq,
                               'doc': trinity.doc})
        ntranscript += 1
        print('n:', ntranscript)

    print('{} transcripts loaded'.format(ntranscript))

    return


def transcriptLoad2(kdb, indent=2):
    """-----------------------------------------------------------------------------------------------------------------
                id TEXT,
                component INTEGER,
                gene INTEGER,
                isoform INTEGER,
                seq TEXT,
                doc TEXT

    -----------------------------------------------------------------------------------------------------------------"""
    transcript_file = input('\n{0}{1}: '.format(' ' * indent, 'Transcript file to load'))
    trinity = Trinity()
    trinity.open(transcript_file)
    ntranscript = 0
    sql = '''
        INSERT INTO transcript
        VALUES ( :id, :component, :gene, :isoform, :seq, :doc )
        '''
    kdb.db.execute("PRAGMA synchronous = OFF")
    while trinity.next():
        kdb.db.execute(sql, {'id': trinity.id, 'component': trinity.component, 'gene': trinity.gene,
                             'isoform': trinity.isoform, 'seq': trinity.seq, 'doc': trinity.doc})
        ntranscript += 1
        print('n:', ntranscript)

    print('{} transcripts loaded'.format(ntranscript))

    return


def transcriptLoadChunk(kdb, indent=2, chunk=1000):
    """
    Load transcripts from file in blocks using loadFromList()
    :param kdb: kollemadb object
    :param indent: number of spaces to indent in output
    :param chunk: number of row to load at once
    :return:
    """
    transcript_file = input('\n{0}{1}: '.format(' ' * indent, 'Transcript file to load'))
    trinity = Trinity()
    trinity.open(transcript_file)

    columns = ['id', 'component', 'gene', 'isoform', 'seq', 'doc']
    ntranscript = 0
    data = []
    kdb.db.execute("PRAGMA synchronous = OFF")
    while trinity.next():
        ntranscript += 1
        if ntranscript % chunk:
            data.append((trinity.id, trinity.component, trinity.gene, trinity.isoform, trinity.seq, trinity.doc))
        else:
            kdb.loadFromList('transcript', columns, data)
            data = []

        print('n:', ntranscript)

    print('{} transcripts loaded'.format(ntranscript))


"""---------------------------------------------------------------------------------------------------------------------
kollema main program
---------------------------------------------------------------------------------------------------------------------"""
version = '0.0.1'
print('Welcome to Kollema, v{0}'.format(version))

root = os.getcwd() + '/.kollema/'
dbfile = root + 'kollema.sql'
if not os.path.isdir(root):
    os.makedirs(root)

kdb = Kollemadb(dbfile=dbfile, new=True)

menu = Menu()
menu.add('main', {'S': 'Select project', 'N': 'New project', 'T': 'Tasks', 'L': 'Load transcripts', 'Q': 'Quit'})
menu.addDispatch('main', {'S': kdb.get, 'N': kdb.fromTerm, 'L': transcriptLoadChunk, 'S': projectSelect})
menu.addTitle('main', 'Current Projects')

select = 'init'
project = 'Not selected'
while True:
    # main event loop
    menu.clear()
    menu.status({'Database file': dbfile, 'Project': project})

    nprojects = kdb.get('project')
    if nprojects:
        menu.clear()
        print(menu.title['main'])
        print(kdb.asFormatted())
    else:
        print('No current projects')

    select = menu.ask('main', 2)
    if select == 'Q': break

    # menu.clear()
    if select == 'N':
        response = menu.dispatch()('project')

    elif select == 'S':
        # select a project from existing projects
        project = menu.dispatch()(kdb)

    elif select == 'L':
        if project:
            response = menu.dispatch()(kdb)
        else:
            print('Project must be selected before loading')

    else:
        print('Unknown option ({})'.format(select))

print('\nThank you')
exit(0)

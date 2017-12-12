"""---------------------------------------------------------------------------------------------------------------------
kollema
Transcriptome annotation
---------------------------------------------------------------------------------------------------------------------"""
import os
import sqlite3 as sql
from kollemadb import Kollemadb
from menu import Menu
from trinity import Trinity


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


def transcriptLoad():
    """-----------------------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------------------"""
    pass
    return


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
menu.addDispatch('main', {'S': kdb.get, 'N': kdb.fromTerm, 'L': transcriptLoad, 'S': projectSelect})
menu.addTitle('main', 'Current Projects')

select = 'init'
project = ''
while True:
    # main event loop

    nprojects = kdb.get('project')
    if nprojects:
        menu.clear()
        print(menu.title['main'])
        print(kdb.asFormatted())
    else:
        print('No current projects')

    select = menu.ask('main', 2)
    if select == 'Q': break

    #menu.clear()
    if select == 'N':
        response = menu.dispatch()('project')

    elif select == 'S':
        # select a project from existing projects
        project = menu.dispatch()(kdb)

    elif select == 'L':
        if project:
            response = menu.dispatch()()
        else:
            print('Project must be selected before loading')

    else:
        print('Unknown option ({})'.format(select))

print('\nThank you')
exit(0)

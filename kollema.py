'''--------------------------------------------------------------------------------
kollema
Transcriptome annotation
--------------------------------------------------------------------------------'''
import os
import sqlite3 as sql
from kollemadb import Kollemadb
from menu import Menu
from trinity import Trinity


def status():
    print('status activated')
    pass


version = '0.0.1'
print('Welcome to Kollema, v{0}'.format(version))

root = os.getcwd() + '/.kollema/'
dbfile = root + 'kollema.sql'
if not os.path.isdir(root):
    os.makedirs(root)

kdb = Kollemadb(dbfile=dbfile, new=True)

menu = Menu()
menu.add('main', {'S': 'Select project', 'N': 'New project', 'T': 'Tasks', 'L': 'Load transcripts', 'Q': 'Quit'})
menu.addDispatch('main', {'S': kdb.get, 'N': kdb.fromTerm})
menu.addTitle('main', 'Current Projects')

select = 'init'
while True:

    nprojects = kdb.get('project')
    if nprojects:
        menu.clear()
        print(menu.title['main'])
        print(kdb.asFormatted())
    else:
        print('No current projects')

    select = menu.ask('main', 2)
    if select == 'Q': break

    menu.clear()
    response = menu.dispatch()('project')

print('\nThank you')
exit(0)

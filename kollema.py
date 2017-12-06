'''--------------------------------------------------------------------------------
kollema
Transcriptome annotation
--------------------------------------------------------------------------------'''
import sqlite3 as sql
from kollemadb import Kollemadb
from menu import Menu


def status():
    print('status activated')
    pass


version = '0.0.1'
print('Welcome to Kollema, v{0}'.format(version))

dbfile = 'kollema.sql'
kdb = Kollemadb(new=True)

menu = Menu()
menu.add('top', {'S': 'Status', 'N':'New project', 'T': 'Tasks', 'L': 'Load transcripts', 'Q': 'Quit'})
x = 1
menu.addDispatch('top', {'S': kdb.get, 'N':kdb.fromTerm})

kdb = Kollemadb()
while True:
    select = menu.ask('top', 2)
    response = menu.dispatch()('project')
    print('response:', response)
    if select == 'Q': break

print('\nThank you')
exit(0)

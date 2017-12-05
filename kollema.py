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


dbfile = 'kollema.sql'
version = '0.0.1'
print('Welcome to Kollema, v{0}'.format(version))

menu = Menu()
menu.add('top',{'S':'Status', 'T':'Tasks', 'L':'Load transcripts', 'Q':'Quit'})

kdb = Kollemadb()
while True:
    select = menu.ask('top', 2)
    if select == 'Q': break

print('\nThank you')
exit(0)

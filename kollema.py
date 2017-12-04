'''--------------------------------------------------------------------------------
kollema
Transcriptome annotation
--------------------------------------------------------------------------------'''
import sqlite3 as sql
from kollemadb import Kollemadb


def menu(id, indent=2):
    '''
    print a menu string using input, and return the response
    :param id: name of the menu
    :param indent: number of spaces to indent, default=2
    :return: selected response
    '''
    response = input('\n{0}{1}: '.format(' ' * indent, menudefs[id]))
    return response.upper()[0]

    # return response[0]

def status():
    print('status activated')
    pass

menudefs = {'top': 'S)tatus  T)asks  L)oad transcripts  Q)uit'}
dispatch = {'S':status()}

dbfile = 'kollema.sql'
version = '0.0.1'
print('Welcome to Kollema, v{0}'.format(version))

kdb = Kollemadb()
while True:
    select = menu('top', 2)
    if select == 'Q': break
    print('  {0} selected'.format(select))

    try:
        dispatch[select]

    except KeyError:
         print('Unknown function selected ({0}'.format(select))

exit(0)

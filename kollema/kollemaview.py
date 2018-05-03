"""=================================================================================================


Michael Gribskov     02 May 2018
================================================================================================="""
import os
import sqlite3 as sq3

import cherrypy
from cherrypy.lib.static import serve_file


class KollemaCherry:
    """=============================================================================================
    Kollema annotation information server

    ============================================================================================="""

    @cherrypy.expose
    def index(self):
        print('path', path)
        return serve_file(os.path.join(static, 'kollema.html'))

    @cherrypy.expose
    def signon(self, new):
        return serve_file(os.path.join(static, 'signon.html'))

    """
    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getData(self, dbfile='gff.db'):
        dbh = sq3.connect(dbfile)
        db = dbh.cursor()

        sql = 'SELECT DISTINCT feature FROM gff;'
        db.row_factory = sq3.Row
        db.execute(sql)
        featurelist = []
        for row in db:
            for key in row.keys():
                if row['feature']:
                    featurelist.append(row['feature'])

        features = '{}'.format(','.join(featurelist))

        return {
            'foo': featurelist
        }

    @cherrypy.expose
    def test(self, type, dbfile='gff.db'):
        dbh = sq3.connect(dbfile)
        db = dbh.cursor()
        sql = 'SELECT attribute FROM gff WHERE feature="{}";'.format(type)
        db.row_factory = sq3.Row
        db.execute(sql)

        idlist = []
        for row in db:
            if row['attribute']:
                field = row['attribute'].split(';')
                for f in field:
                    if f.startswith('ID='):
                        f = f.replace('ID=','')
                        if f not in idlist:
                            idlist.append(f)
                        break

        response = ''
        for id in idlist:
            response += '<option value="{}">{}</option>\n'.format(id,id)

        return response
    """


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # path to root directory
    path = os.path.abspath(os.path.dirname(__file__))
    static = os.path.join(path, 'static')
    css = os.path.abspath('./static/css')
    print('css:', css)
    images = os.path.join(static, 'images')
    fav = os.path.join(static, 'images', 'kappa-64.png')
    print('fav:', fav)

    config = {
        'global': {
            'server.socket_host': '127.0.0.1',
            'server.socket_port': 8080,
            'server.thread_pool': 2
        },
        '/css':{
            'tools.staticdir.on': True,
            'tools.staticdir.dir': css
        },
        '/images': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': images
        },
        '/favicon.png': {
            'tools.staticfile.on': True,
            'tools.staticfile.filename': fav
        }

    }

    cherrypy.quickstart(KollemaCherry(), '/', config)

    cherrypy.engine.exit()

    exit(0)

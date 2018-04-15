"""=================================================================================================
cherryPy app with ajax/jquery

Michael Gribskov     15 April 2018
================================================================================================="""
import os
import sqlite3 as sq3

import cherrypy
from cherrypy.lib.static import serve_file


class GffCherry:
    """=============================================================================================
    Gff information server

    ============================================================================================="""

    @cherrypy.expose
    def index(self):
        return serve_file(os.path.join(path, 'index.html'))

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getData(self, dbfile='gff.db'):
        dbh = sq3.connect(dbfile)
        db = dbh.cursor()

        sql = 'SELECT DISTINCT feature FROM gff;'
        db.row_factory = sq3.Row
        db.execute(sql)
        featurelist = []
        features = ''
        for row in db:
            for key in row.keys():
                if row['feature']:
                    featurelist.append(row['feature'])

        features = ','.join(featurelist)

        return {
            'foo': features,
            'baz': 'another one'
        }


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    path = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(path, 'html')

    config = {
        'global': {
            'server.socket_host': '127.0.0.1',
            'server.socket_port': 8081,
            'server.thread_pool': 8
        }
    }

    cherrypy.quickstart(GffCherry(), '/', config)


    cherrypy.engine.exit()

    exit(0)

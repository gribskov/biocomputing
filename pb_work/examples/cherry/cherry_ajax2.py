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
        return serve_file(os.path.join(path, 'index2.html'))

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

        return {
            'foo': featurelist,
            'baz': 'another one' }

    @cherrypy.expose
    # @cherrypy.tools.json_in()
    # @cherrypy.tools.json_out()
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
                        f = f.replace('ID=', '')
                        if f not in idlist:
                            idlist.append(f)
                        break

        response = ''
        for id in idlist:
            response += '<option value="{}">{}</option>\n'.format(id, id)

        print(type, response)
        return response


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    path = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(path, 'html')

    config = {
        'global': {
            'server.socket_host': '127.0.0.1',
            'server.socket_port': 8080,
            'server.thread_pool': 8
        }
    }

    cherrypy.quickstart(GffCherry(), '/', config)

    cherrypy.engine.exit()

    exit(0)

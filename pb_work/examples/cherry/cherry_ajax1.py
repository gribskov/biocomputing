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
        return serve_file(os.path.join(path, 'index0.html'))

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getData(self):
        return {
            'foo': 'bar sent from cherryPy',
            'baz': 'faz - another one'
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
            'server.socket_port': 8080,
            'server.thread_pool': 8
        }
    }

    cherrypy.quickstart(GffCherry(), '/', config)

    cherrypy.engine.exit()

    exit(0)

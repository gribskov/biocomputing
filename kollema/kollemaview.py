"""=================================================================================================
Kollema annotation system

Michael Gribskov     02 May 2018
================================================================================================="""
import os
import sqlite3 as sq3
from kollemadb import Kollemadb
from trinity.trinity import Trinity

import cherrypy
from cherrypy.lib.static import serve_file


class KollemaCherry:
    """=============================================================================================
    Kollema annotation information server

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        constuctor
        -----------------------------------------------------------------------------------------"""
        self.dbfile = "kollema.sqlite3"
        self.user = 'test'

    @cherrypy.expose
    def index(self):
        """-----------------------------------------------------------------------------------------
        Splash page
        :return: kollema.htm, static page
        -----------------------------------------------------------------------------------------"""
        # print('path', path)
        return serve_file(os.path.join(static, 'kollema.html'))

    @cherrypy.expose
    def dashboard(self, firstname=None, lastname=None, phone=None, email=None):
        """-----------------------------------------------------------------------------------------
        The dashboard controls most display and analyses through ajax interactions.
        :param firstname:
        :param lastname:
        :param phone:
        :param email:
        :return: dashboard.html, static page
        -----------------------------------------------------------------------------------------"""
        kdb = Kollemadb(self.dbfile)

        result = kdb.getByUser('user', email)
        if len(result) == 0:
            # no matching user in user table: unknown user or new
            print('row is none')
            if firstname == '':
                # unknown user
                # TODO error handling when data is missing, javascript?
                pass
            else:
                # new user
                kdb.set('user', {'firstname': firstname, 'lastname': lastname, 'phone': phone,
                                 'email': email})
                self.user = email
                cherrypy.session['user'] = email
        else:
            # match in user table: authenticated
            for row in result:
                for key in row.keys():
                    print(row[key], end='\t')
                print()

            self.user = email
            cherrypy.session['user'] = email

        # create html page, substitute the session user for the token $$user in the static html
        dashboard = open('static/dashboard.html', 'r').read()

        return dashboard.replace('$$user', cherrypy.session['user'])

    # end of dashboard

    @cherrypy.expose
    def getProjects(self, user):
        """-----------------------------------------------------------------------------------------
        Ajax function for project list
        :return:
        -----------------------------------------------------------------------------------------"""
        kdb = Kollemadb(self.dbfile)
        result = kdb.getByUser('projects', user)

        html = '<h6>Projects: user {}</h6><br>\n'.format(user)
        html += '<table id="project_table" style="font-size:1.0rem;">\n'
        html += '<tr><th>Name</th><th>Description</th><th>Created</th><th>Status</th><th>Updated</th>\n'
        if len(result) == 0:
            html += '<tr><td colspan="5">No projects found for {}</td></tr>\n'.format(user)

        else:
            for row in result:
                html += '<tr><td class="sm-medium">{}</td><td class="sm-wide">{}</td><td class="sm-narrow">{}</td><td class="sm-narrow">{}</td><td class="sm-narrow">{}</td></tr>\n'.format(
                    row['name'], row['description'], row['created'], row['status'], row['updated'])

        html += '</table>'
        # print('getProjects-html: ', html)

        return html

    # end of getprojects

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def addProject(self, name=None, desc=None):
        """-----------------------------------------------------------------------------------------
        Add a new project to the project table

        :param name: string, project name
        :param desc: string, project description
        :return: JSON object, user_id
        -----------------------------------------------------------------------------------------"""
        kdb = Kollemadb(self.dbfile)

        # TODO escape newline so it is stored in db
        desc.replace('\n', '\\n')
        kdb.set('projects', {'email': self.user, 'name': name, 'description': desc})

        return {'user': self.user}

    # end of addProject

    @cherrypy.expose
    def setProject(self, project_id=None):
        """-----------------------------------------------------------------------------------------
        look up selected project and populate screens
        -----------------------------------------------------------------------------------------"""
        kdb = Kollemadb(self.dbfile)
        cherrypy.session['project_id'] = project_id

        return "Under development:{}".format(project_id)

    @cherrypy.expose
    def transcriptLoad(self, transcript_file=None, override=None):
        """-----------------------------------------------------------------------------------------
         CREATE TABLE `transcript` (
            `transcript_id`	INTEGER PRIMARY KEY AUTOINCREMENT,
            `project_id`	INTEGER,
            `email`	TEXT,
            `cluster`	INTEGER,
            `component`	INTEGER,
            `gene`	INTEGER,
            `isoform`	INTEGER,
            `seq`   TEXT,
            `doc`   TEXT,
            `path`	TEXT,
            `timestamp`	TEXT
        );
        -----------------------------------------------------------------------------------------"""
        try:
            project_id = cherrypy.session['project_id']
        except KeyError:
            project_id = None

        if project_id is None:
            print('Kollemaview.transcriptLoad - project is undefined')
            return
        user = cherrypy.session['user']

        kdb = Kollemadb(self.dbfile)
        trinity = Trinity()
        trinity.fh = transcript_file.file
        sql = '''
            INSERT INTO transcript
            (project_id, email, cluster, component, gene, isoform, seq, doc)
            VALUES ( :project_id, :user_id, :cluster, :component, 
                     :gene, :isoform, :seq, :doc )
            '''
        kdb.db.execute('PRAGMA synchronous = NORMAL')
        kdb.db.execute('PRAGMA journal_mode=WAL')
        kdb.dbh.set_trace_callback(print)

        ntranscript = 0
        while trinity.next():
            cluster, component, gene, isoform = Trinity.splitID(trinity.id)
            kdb.db.execute(sql, {'project_id': project_id,
                                 'user_id': user,
                                 'cluster': cluster,
                                 'component': int(component),
                                 'gene': int(gene),
                                 'isoform': int(isoform),
                                 'seq': trinity.seq,
                                 'doc': trinity.doc
                                 })
            ntranscript += 1

        print('Kollemaview.transcriptLoad - {} transcripts loaded'.format(ntranscript))

        return


# --------------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    kollema = KollemaCherry()

    # path to root directory
    path = os.path.abspath(os.path.dirname(__file__))
    static = os.path.join(path, 'static')
    css = os.path.abspath('./static/css')
    # print('css:', css)
    images = os.path.join(static, 'images')
    fav = os.path.join(static, 'images', 'kappa-64.png')
    # print('fav:', fav)

    config = {
        'global': {
            'server.socket_host': '127.0.0.1',
            'server.socket_port': 8080,
            'server.thread_pool': 4
        },
        '/static': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': static

        },
        '/css': {
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
        },
        '/': {
            'tools.sessions.on': True
        }

    }

    cherrypy.quickstart(KollemaCherry(), '/', config)

    cherrypy.engine.exit()

    exit(0)

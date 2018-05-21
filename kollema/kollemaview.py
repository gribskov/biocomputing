"""=================================================================================================
Kollema annotation system

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
        dbh = sq3.connect(self.dbfile)
        db = dbh.cursor()

        # print('first:{}\tlast:{}\temail:{}\tphone:{}'.format(firstname, lastname, email, phone))

        # known user
        sql = 'SELECT * FROM user WHERE email="{}"'.format(email)
        # print(sql)
        db.row_factory = sq3.Row
        db.execute(sql)
        row = db.fetchone()
        if row is None:
            # unknown user or new
            if firstname == '':
                # unknown
                pass
            else:
                # new user
                sql = 'INSERT INTO user VALUES ( {}, "{}", "{}", "{}", "{}", {});'.format(
                    'Null', firstname, lastname, phone, email, 'Null')
                # print(sql)
                db.execute(sql)
                dbh.commit()
                self.user = email
                cherrypy.session['user'] = email

        else:
            # authenticated
            for row in db:
                for key in row.keys():
                    print(row[key], end='\t')
                print()

            self.user = email
            cherrypy.session['user'] = email

        # create html page
        dashboard = open('static/dashboard.html', 'r').read()

        return dashboard.replace('$$user', cherrypy.session['user'])

    # end of dashboard

    @cherrypy.expose
    def getProjects(self, user):
        """-----------------------------------------------------------------------------------------
        Ajax function for project list
        :return:
        -----------------------------------------------------------------------------------------"""
        dbh = sq3.connect(self.dbfile)
        db = dbh.cursor()
        sql = 'SELECT * FROM projects WHERE email="{}"'.format(user)
        db.row_factory = sq3.Row
        db.execute(sql)
        result = db.fetchall()

        html = '<span id="utitle">User: {}</span><br>\n'.format(user)
        html += '<table id="project_table">\n'
        # html += '<col class="sm-medium">\n'
        # html += '<col class ="sm-wide">\n'
        # html += '<col class ="sm-narrow">\n'
        # html += '<col class ="sm-narrow">\n'
        # html += '<col class ="sm-narrow">\n'
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
        :return:
        -----------------------------------------------------------------------------------------"""
        dbh = sq3.connect(self.dbfile)
        db = dbh.cursor()
        # TODO escape newline so it is stored in db
        sql = 'INSERT INTO projects VALUES ( {}, "{}", "{}", "{}", "{}", "{}", "{}", "{}" )'.format(
            'Null', self.user, name, desc, '', '', '', '')
        print(sql)

        db.row_factory = sq3.Row
        db.execute(sql)
        dbh.commit()

        # print('addProject-user:', self.user )
        return { 'user': self.user }

    # end of addProject


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
            'server.socket_port': 8081,
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

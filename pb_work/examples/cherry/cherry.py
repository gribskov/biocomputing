import time
import cherrypy
import sqlite3 as sq3


class seqServe():
    """=============================================================================================
    cherrypy sequence analysis
    ============================================================================================="""

    @cherrypy.expose
    def test_db(self, dbfile='gff.db'):
        """-----------------------------------------------------------------------------------------
        check to see the database works
        :return:
        -----------------------------------------------------------------------------------------"""
        dbh = sq3.connect(dbfile)
        db = dbh.cursor()

        gff = self.page_header()
        gff += '<pre>\n'

        sql = 'SELECT name FROM my_db.sqlite_master WHERE type="table";'
        db.row_factory = sq3.Row
        db.execute(sql)
        for row in db:
            for key in row.keys():
                print('{}:{}'.format(key, row[key]), end='\t')
            print()

        gff += '</pre>\n'
        gff += self.page_footer()

        return gff

    @cherrypy.expose
    def load(self):
        """-----------------------------------------------------------------------------------------
        Query page to get and load a gff file
        :return: html page via cherrypy
        -----------------------------------------------------------------------------------------"""
        load = self.page_header()
        load += self.file_form('gffread')
        load += self.page_footer()

        return load

    @cherrypy.expose
    def gffread(self, gff_file, dbfile='gff.db', samplesize=10):
        """-----------------------------------------------------------------------------------------
        display 10 lines of the uploaded gff file
        :param gff_file:
        :return:
        -----------------------------------------------------------------------------------------"""
        dbh = sq3.connect(dbfile)
        db = dbh.cursor()

        gff = self.page_header() + self.title()
        gff += '<br/>Sample of the first {} features\n<br/>\n'.format(samplesize)
        gff += '<pre>\n'

        timestart = time.time()
        nline = 0
        for line in gff_file.file:
            line = line.decode()
            if line.startswith('#'):
                continue

            field = line.split('\t')

            # fix data types
            field[3] = int(field[3])
            field[4] = int(field[4])
            if field[5] == '.':
                field[5] = 0.0
            else:
                field[5] = float(field[5])

            sql = 'INSERT INTO gff VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ? )'
            db.execute(sql, tuple(field[:]))

            if nline < samplesize:
                gff += line
            nline += 1

        dbh.commit()
        timeend = time.time()

        gff += '</pre>'
        gff += '{} features uploaded in {:.1f} seconds\n'.format(nline, timeend - timestart)
        gff += self.page_footer()

        return gff

    @cherrypy.expose
    def query(self):
        """-----------------------------------------------------------------------------------------
        Query page for genes
        :return:
        -----------------------------------------------------------------------------------------"""
        query = self.page_header()
        query += '<h2><hr/>Query<hr/></h2>\n'
        query += '    <form action="draw" method="post"><br/>\n'
        query += '        <input type="text" name="gene" value="AT4G00090"><br/>\n'
        query += '        <input type="submit" value="Search">\n'
        query += '    </form><br/>\n'

        query += self.page_footer()
        return query

    @cherrypy.expose
    def draw(self, gene, dbfile='gff.db'):
        """-----------------------------------------------------------------------------------------
        make sql gene query and draw gene
        :return:
        -----------------------------------------------------------------------------------------"""
        dbh = sq3.connect(dbfile)
        db = dbh.cursor()

        draw = self.page_header()
        draw += '<h2><hr/>Gene {}<hr/></h2>\n'.format(gene)

        sql = 'SELECT * FROM gff WHERE attribute LIKE "%{}%" AND feature == "CDS";'.format(gene)
        print('sql:',sql)
        select = db.execute(sql)

        draw += '    <pre>\n'
        for row in select:
            draw += str(row) + '\n'

        draw += '    </pre>\n'
        draw += self.page_footer()

        return draw

        return draw

    def page_header(self):
        """-----------------------------------------------------------------------------------------
        HTML header.  open html and body, includes complete head element
        :return:
        -----------------------------------------------------------------------------------------"""
        return '<html>\n<head>\n</head>\n<body>\n<br/>' + self.menu_main() + self.title()

    def page_footer(self):
        """-----------------------------------------------------------------------------------------
        HTML footer, closes boday and html
        :return:
        -----------------------------------------------------------------------------------------"""
        return '\n</body\n</html\n'

    def menu_main(self):
        """-----------------------------------------------------------------------------------------
        main navigation menu
        :return:
        -----------------------------------------------------------------------------------------"""
        style = 'display: inline;\
padding: 0.7em 1em 0.7em 1em;\
border: 1px solid #FFFFFF;\
border-radius: 5px;\
margin:10px 0px 10px 0px;\
background-color: #000088;\
color: #FFFFFF;text-decoration:none;'

        menu = '\n<!- main menu - seqServe.menu_main ->\n'
        menu += '    <a href="load" style="{}">load</a>\n'.format(style)
        menu += '    <a href="query" style="{}">query</a>\n'.format(style)
        return menu

    def title(self, margin=4, indent=4):
        """-----------------------------------------------------------------------------------------
        Standard HTML title for all pages
        :return:
        -----------------------------------------------------------------------------------------"""
        level = 0
        space = ' ' * (margin + level * indent)
        return '{}<h1>MyGene</h1><br/>\n'.format(space)

    def file_form(self, action, margin=4, indent=4):
        """-----------------------------------------------------------------------------------------
        html to get a file as multipart/form-data
        :return:
        -----------------------------------------------------------------------------------------"""
        level = 0
        space = ' ' * (margin + level * indent)
        html = '{}<form action="{}" method="post" enctype="multipart/form-data">\n'.format(space,
                                                                                           action)

        level += 1
        space = ' ' * (margin + level * indent)
        html += '{}Enter filename <input type="file" name="gff_file"/><br/>\n\n'.format(space)
        html += '{}<input type="submit" value="Upload">\n'.format(space)

        level -= 1
        space = ' ' * (margin + level * indent)
        html += '{}</form>\n'.format(space)

        return html


# ==================================================================================================
# main/testing
# ==================================================================================================
if __name__ == '__main__':
    cherrypy.config.update({'server.socket_port': 8082, })
    cherrypy.quickstart(seqServe())

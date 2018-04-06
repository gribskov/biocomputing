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

        gff = self.page_header() + self.title()
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
        load = self.page_header() + self.title()
        load += self.file_form('gffread')
        load += self.page_footer()

        return load

    @cherrypy.expose
    def gffread(self, gff_file, dbfile='gff.db'):
        """-----------------------------------------------------------------------------------------
        display 10 lines of the uploaded gff file
        :param gff_file:
        :return:
        -----------------------------------------------------------------------------------------"""
        dbh = sq3.connect(dbfile)
        db = dbh.cursor()


        gff = self.page_header() + self.title()
        gff += '<pre>\n'
        nline = 0
        for line in gff_file.file:
            line = line.decode()
            if line.startswith('#'):
                continue

            gff += line
            nline += 1

            field = line.split('\t')

            # fix data types
            field[3] = int(field[3])
            field[4] = int(field[4])
            if field[5] == '.':
                field[5] = 0.0
            else:
                field = float(field[5])


            sql = 'INSERT INTO gff VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ? )'
            db.execute(sql, tuple(field[:]))

            if nline > 9:
                break

        dbh.commit()
        gff += '</pre>'
        gff += self.page_footer()

        return gff

    def page_header(self):
        """-----------------------------------------------------------------------------------------
        HTML header.  open html and body, includes complete head element
        :return:
        -----------------------------------------------------------------------------------------"""
        return '<html>\n<head>\n</head>\n<body>\n'

    def page_footer(self):
        """-----------------------------------------------------------------------------------------
        HTML footer, closes boday and html
        :return:
        -----------------------------------------------------------------------------------------"""
        return '\n</body\n</html\n'

    def title(self, margin=4, indent=4):
        """-----------------------------------------------------------------------------------------
        Standard HTML title for all pages
        :return:
        -----------------------------------------------------------------------------------------"""
        level = 0
        space = ' ' * (margin + level * indent)
        return '{}<h1>MyGene</h1>\n'.format(space)

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
    cherrypy.config.update({'server.socket_port': 8081, })
    cherrypy.quickstart(seqServe())

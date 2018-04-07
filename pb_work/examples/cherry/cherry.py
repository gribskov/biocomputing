import time
import os
import cherrypy
import sqlite3 as sq3
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Rectangle


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

        sql = 'SELECT * FROM gff WHERE attribute LIKE "%ID={}%" AND feature == "CDS";'.format(gene)
        print('sql:', sql)
        select = db.execute(sql)

        draw += '    <pre>\n'
        exonlist = []
        for row in select:
            exonlist.append(row)
            draw += str(row) + '\n'

        img = self.gene_img(exonlist)
        draw += '\n<img src="{}">\n'.format(img)

        draw += '    </pre>\n'
        draw += self.page_footer()

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

    def gene_img(self, elements):
        """-----------------------------------------------------------------------------------------
        Draw the gene exons/introns as boxes/lines using matplotlib
        :param elements: list of tuples
        :return: image file name
        -----------------------------------------------------------------------------------------"""
        majorlocator = MultipleLocator(1000)
        minorlocator = MultipleLocator(100)
        majorformatter = FormatStrFormatter('%d')

        minpos = 1e8
        maxpos = 0
        # find minimum and maximum positions
        for feature in elements:
            minpos = min(minpos, int(feature[4]))
            maxpos = max(maxpos, int(feature[5]))

        begin = minpos - 500
        end = maxpos + 500

        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_subplot(111)
        ticks = [i for i in range(begin, end) if i % 100 == 0]

        plt.xticks(ticks)
        plt.xlim(begin, end)
        plt.ylim(0, 100)

        plt.plot([minpos, maxpos], [50.0, 50.0], color='black', linewidth=2.0)
        y1 = 45.0
        y2 = 55.0
        for feature in elements:
            x1 = float(feature[4])
            x2 = float(feature[5])
            plt.fill([x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], color='red', zorder=3)
            plt.plot([x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], color='black', linewidth=1.5, zorder=4)

        ax.xaxis.set_major_locator(majorlocator)
        ax.xaxis.set_major_formatter(majorformatter)
        ax.xaxis.set_minor_locator(minorlocator)

        fig.savefig('img/plot.png')
        return 'img/plot.png'


# ==================================================================================================
# main/testing
# ==================================================================================================
if __name__ == '__main__':
    file_path = os.getcwd()

    cherrypy.server.socket_host = "127.0.0.1"
    cherrypy.server.socket_port = 8080

    cherrypy.quickstart(seqServe(), '/', {
        "/img": {
            "tools.staticdir.on": True,
            "tools.staticdir.dir": os.path.join(file_path, "img"),
        }
    })

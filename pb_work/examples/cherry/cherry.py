import cherrypy


class seqServe():
    """=============================================================================================
    cherrypy sequence analysis
    ============================================================================================="""

    @cherrypy.expose
    def home(self):
        return 'welcome home'

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
    def gffread(self, gff_file):
        """-----------------------------------------------------------------------------------------
        display 10 lines of the uploaded gff file
        :param gff_file:
        :return:
        -----------------------------------------------------------------------------------------"""
        x = gff_file

        gff = self.page_header() + self.title()
        gff += '<pre>\n'
        nline = 0
        for line in gff_file.file:
            gff += line.decode()
            nline += 1
            if nline > 9:
                break
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
    cherrypy.quickstart(seqServe())

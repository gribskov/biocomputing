import re
import sys


class Blast(object):
    """=============================================================================================
    Blast class is intended to iterate over a file of blast results.Can also be 
    used for other programs with blast-like output such as diamond.
    
    Multiple reading functions can be written and selected as the active reading
    function using blast.read( function ).  Currently, only a tabular reading
    function intended for use with the output from Diamond is implemented.

    execute the class from the command line for testing

    usage

    blast = Blast()
    blast.open( 'file.blastx' )
    while blast.next():
        blast.dump()        # print out current content of object
        ... do something

    or blast.read()     # read a line using the selected internal read function

    ============================================================================================="""

    def __init__(self, file=''):
        """-----------------------------------------------------------------------------------------
        initialize data structure.  The following fields should always be defined:
            file, fh    # filename, filehandle
            read        # read function for use in next
            format      # space delimited string of fields in tabular output
            fields      # ordered list of fields in tabular output, automatically 
                        # created from format. Default is diamond outfmt=6
            line        # last string read
        In addition, variable attributes depending on the format will be filled.  
        In the case of tabulara format, the attributes will have the names given i
        in the fields list.
        usage
            blast = Blast()
        -----------------------------------------------------------------------------------------"""
        self.file = ''
        self.fh = None
        self.line = ''
        self.read = self.readTabular
        self.fields = []
        self.format = \
            'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

        # if file is provided, open it
        if file:
            self.new(file)

    def dump(self, indent=4):
        """-----------------------------------------------------------------------------------------
        print the attributes of the object
        usage
            blast.dump()
            blast.dump(5)
            blast.dump(indent=5)
        -----------------------------------------------------------------------------------------"""
        classname = re.search("__main__\.([^']+)'>", str(self.__class__))
        classname = classname.group(1)
        print('{0}dumping{1}'.format(' ' * indent, classname))
        for attr in dir(self):
            if hasattr(self, attr):
                if attr[0] == '_':
                    continue
                descr = str(getattr(self, attr))
                if descr[0] == '<':
                    continue
                print('{0}{1}.{2} = {3}'.format(' ' * indent, classname, attr, descr))

    def next(self):
        """-----------------------------------------------------------------------------------------
        retrieve the next entry from the file, using the currently registered read method
        usage
            while blast.read:
               ... do something
        -----------------------------------------------------------------------------------------"""
        return self.read()

    def new(self, file):
        """-----------------------------------------------------------------------------------------
        open a result file for reading and store the filehandle
        usage
        blast.new('file.blastx')
        -----------------------------------------------------------------------------------------"""
        try:
            fh = open(file, 'r')

        except:
            print('Blast: could not open file ({0})'.format(file), file=sys.stderr)
            return False

        self.file = file
        self.fh = fh
        return True

    def readTabular(self):
        """-----------------------------------------------------------------------------------------
        Read one hit in tabular format.  Diamond does not provide any column ID information
        so it must be provided via blast.format(). Generally this function should be called
        via the blast.read method.
        Returns the line read from the file
        usage
            blast.readTabular()
        -----------------------------------------------------------------------------------------"""
        line = '#'
        while line.startswith('#'):
            # TODO needs to be improved to record the names of sequences with no hits
            line = self.fh.readline()

        if line:
            line = line.rstrip()
            # print('line:', line)
            token = line.split()
            # TODO check to see fields are defined
            n = 0
            last_key = ''
            for key in self.fields:
                self.__dict__[key] = token[n]
                last_key = key
                # print('{0} => {1}'.format(key,token[n]))
                n += 1

            if n < len(token):
                self.__dict__[last_key] += ' ' + ' '.join(token[n:])

            self.line = line
            return line
        else:
            # end of file
            return False

    def setFormat(self, fmt='', preset='diamond_doc'):
        """-----------------------------------------------------------------------------------------
        change the fields in the read format. fields is a list of the attributes available for this
        Blast search.  A instance attribute is created for each known attribute.

        diamond_doc format = 'qid qlen qstart qend sid slen start send allen pid score evalue doc'
                             for adding annotation to transcripts

        usage
            nfields = blast.setFormat('qid sid qcov pid len evalue')
            nfields = blast.setFormat(preset='diamond_doc')
        -----------------------------------------------------------------------------------------"""
        # first delete all the existing fields
        for attribute in self.fields:
            del self.__dict__[attribute]

        # now create the new fields
        if fmt:
            self.format = fmt

        else:
            self.fmt = preset

        self.fields = self.format.split()
        for attribute in self.fields:
            self.__dict__[attribute] = None

        return len(self.fields)

    def toDict(self):
        """-----------------------------------------------------------------------------------------
        Return the current line as a dictionary
        :return:
        -----------------------------------------------------------------------------------------"""
        info = {}
        for key in self.fields:
            info[key] = self.__dict__[key]

        return info


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    print('Blast object', end='\n\n')
    blast = Blast()
    print('Blast - test opening non-existant file')
    if blast.new('blast/diamond.blastx'):
        print('    success')
    else:
        print('    failure')
    print('Blast - test opening existing file')
    if blast.new('data/diamond.blastx'):
        print('    success')
    else:
        print('    failure')

    print('')
    print('Blast - format for tabular file')
    print('    default read format:', blast.format)
    nfields = blast.setFormat('qid sid qcov pid len evalue')
    print('    modified read format:', blast.format)
    print('    fields:', nfields, '\tformat:', blast.format)
    for f in blast.fields:
        print('        ', f, '\t', getattr(blast, f))

    print('')
    print('Blast - test reading')
    print('    read with blast.read:', blast.read())
    print('    dump with indent=default (4)')
    blast.dump()
    print('')
    print('    read second line with read:', blast.read())
    print('    dump with indent=5')
    blast.dump(5)

    print('')
    print('Blast - Test next function')
    n = 0
    while blast.next():
        n += 1
        print('   ', n, blast.line)

        # if n > 4: break

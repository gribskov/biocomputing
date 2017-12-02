import re
import sys
from pprint import pprint

class Blast(object):
    '''--------------------------------------------------------------------------
    Blast class is intended to iterate over a file of blast results.Can also be 
    used for other programs with blast-like output such as diamond.

    usage

    blast = Blast()
    blast.open( 'file.blastx' )
    while blast.next():
        ... do something

    --------------------------------------------------------------------------'''

    def __init__(self):
        '''
        initialize data structure.  since the contents will vary, the data 
        structure is a dictionary corresponding to a single result
        usage
            blast = Blast()
        '''
        self.file   = ''
        self.fh     = None
        self.read   = self.readTabular
        self.fields = []
        self.setFormat('qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore')


    def dump(self,indent=4):
       '''
       print the attributes of the object
       usage
           blast.dump()
           blast.dump(5)
           blast.dump(indent=5)
       '''
       myclass = re.search("__main__\.([^']+)'>", str(self.__class__))
       print('dumping', myclass.group(1))
       for attr in dir(self):
           if hasattr( self, attr ):
               if attr[0] == '_': continue
               descr = str(getattr(self, attr))
               if descr[0] == '<': continue
               print( '{0}obj.{1} = {2}'.format(' '*indent, attr, descr))


    def next(self):
        self.read()


    def open(self,file):
        '''
        open a rusult file for reading and store the filehandle
        usage
        blast.open('file.blastx')
        '''
        try:
            fh = open(file,'r')

        except:
           print('Blast: could not open file ({0})'.format(file), file=sys.stderr)
           return False 

        self.file = file
        self.fh   = fh

    def readTabular(self):
        '''
        '''
        line = self.fh.readline()
        #print('line:', line)
        token = line.split()
        n = 0
        for key in self.fields:
            self.__dict__[key] = token[n]
            #print('{0} => {1}'.format(key,token[n]))
            n += 1
        return(line)


    def setFormat(self,fmt):
        '''
        change the fields in the read format
        usage
            nfields = blast.setFormat('qid sid qcov pid len evalue')
        '''
        self.format = fmt
        self.fields = self.format.split()

        return len(self.fields)


if __name__ == '__main__':
    print('Blast object')
    blast = Blast()
    blast.open('blast/diamond.blastx')
    blast.open('data/diamond.blastx')

    print('read format:', blast.format)

    format='qid sid qcov pid len evalue'
    nfields = blast.setFormat(format)
    print('fields:',nfields, 'format:', format)
    print('read:',blast.read())
    blast.dump()
    print('read:',blast.read())
    blast.dump(5)

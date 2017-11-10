'''
Fasta sequence class.  Supports iteration over a multi-fasta file
    id
    documentation
    sequence
'''
class Fasta:
	
	def __init__(self):
		self.filename = ""

	def open(self,filename):
		'''
		open a file for reading
		'''
		self.fh = open(filename, 'r')
		print("opening a file:",filename)

	def next(self);
		'''
		return the next entry from an open file into the object
		'''
		self.read();
		if self.id == ""
			return 0
		else:
			return 1 

	def read(self):
		for line in inp.input():
    		line = line.rstrip('\n')

			print( line )
			if line[0] == '>':
			    fasta.append( { 'id': '', 'documentation': '', 'sequence': '' } )
			    thisfasta = fasta[-1]
			    line = line.lstrip( '>' )
			    try:
			        thisfasta['id'], thisfasta['documentation'] = line.split( " ", 1 )
			    except ValueError:
			        'documentation is missing'
			        thisfasta['id'] = line
			else:
			    thisfasta['sequence'] += line

print( 'fasta' )

# import fileinput as inp

# fasta = []
# for line in inp.input():
#     line = line.rstrip('\n')

#     print( line )
#     if line[0] == '>':
#         fasta.append( { 'id': '', 'documentation': '', 'sequence': '' } )
#         thisfasta = fasta[-1]
#         line = line.lstrip( '>' )
#         try:
#             thisfasta['id'], thisfasta['documentation'] = line.split( " ", 1 )
#         except ValueError:
#             'documentation is missing'
#             thisfasta['id'] = line
#     else:
#         thisfasta['sequence'] += line

# print( fasta )
# print( all )

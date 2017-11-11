'''
Fasta sequence class.  Supports iteration over a multi-fasta file
    id
    documentation
    sequence
'''
buffer = ''

class Fasta:

	def __init__(self):
		self.filename = ""
		self.id = ''
		self.doc = ''
		self.seq = ''

	def open(self,filename):
		'''
		open a file for reading
		'''
		self.fh = open(filename, 'r')
		print("opening a file:",filename)

	def next(self):
		'''
		return the next entry from an open file into the object
		'''
		self.read();
		if not self.id:
			return 0
		else:
			return 1 

	def read(self):
		'''
		read one seqeunce from the file, leave the following line in buffer
		usage:
		fasta.read()
		'''
		if buffer:
			# if buffer has something in it it could be an ID line
			self.id, self.doc = getId(buffer)

		for line in self.fh:
			line = line.rstrip('\n')

			print( line )
			if line[0] == '>':
			    self.id = ""
			    self.doc = ""
			    self.seq = ""
			    line = line.lstrip( '>' )
			    self.id, self.doc = getID(line)
			else:
				self.seq += line

print( 'fasta' )
'''
helper functions, not designed for external use
'''
def getID(line):
	'''
	break a line into ID and documentation and return both
	id will be stripped of >
	documentation will be and empty string if there is nothing following the ID
	'''
	try:
		id, doc = line.split( " ", 1 )		        
	except ValueError:
		'documentation is missing'
		id = line
	return id, doc

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

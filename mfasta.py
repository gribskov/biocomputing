'''
read a multiple fasta file
store in a an array of dictionaries with keys
    id
    documentation
    sequence
'''
import fileinput as inp

fasta = []
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

print( fasta )
print( all )

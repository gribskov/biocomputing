'''
read a single sequence in fasta format from a file. Store information as
id, doc, seq

format
>id documentation
sequence

id is a single token terminated by a space
documentation is any number of tokens terminated by a newline
sequence is any number of lines, terminated by end-of-file
'''
import sys

id = doc = seq = ''
nlines = 0
for line in sys.stdin:
    line = line.rstrip()
    nlines += 1
    if line[0] == '>':
        # a title line
        id,doc = line.split( ' ', maxsplit=1 )
        id = id.lstrip( ' >' )

    else:
        # a sequence line
        seq += line

print( nlines, "lines read" )
print( "id:", id )
print( "doc:", doc )
print( "seq:", seq )


'''
read a single sequence in fasta format from a file

format
>id documentation
sequence

id is a single token terminated by a space
documentation is any number of tokens terminated by a newline
sequence is any number of lines, terminated by end-of-file
'''
file = 'data/1seq.fa'
fasta = open( file, 'r' )

'''most basic method to read files'''
n = 0
for line in fasta:
    n += 1
    print( n, ":", line )

fasta.close

'''preferred method'''
with open( file, 'r' ) as fasta:
    n = 0
    for line in fasta:
        n += 1
        print( n, ":", line, end="" )

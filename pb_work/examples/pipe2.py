"""=================================================================================================
create a pipe between two unix commands using subprocess.Popen
pipe is
grep '>' fasta | wc
which counts the number of sequences in a fasta file
================================================================================================="""
import subprocess as sub

file = '../../data/Trinity.fasta'
fasta = open(file, 'r')

proc = sub.Popen("grep '>'", shell=True, stdin=fasta, stdout=sub.PIPE, stderr=sub.PIPE)
wc = sub.Popen("wc", shell=True, stdin=proc.stdout, stdout=sub.PIPE, stderr=sub.PIPE)

for line in wc.stdout:
    print(line.decode())
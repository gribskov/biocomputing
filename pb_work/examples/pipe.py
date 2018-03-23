"""=================================================================================================
run a simple unix command, grep, on an input file using subprocess.Popen()
================================================================================================="""
import subprocess as sub

file = '../../data/Trinity.fasta'
fasta = open(file, 'r')

proc = sub.Popen("grep '>'", shell=True, stdin=fasta, stdout=sub.PIPE, stderr=sub.PIPE )

for line in proc.stdout:
    print(line.decode())
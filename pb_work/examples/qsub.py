"""=================================================================================================
use subprocess.Popen to submit Torque/PBS job
================================================================================================="""
import subprocess as sub

proc = sub.Popen('qsub', shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE,
                 close_fds=True)

# put the commands in a string
job_string = """#!/bin/bash
#PBS -N testpy
#PBS -q standby
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
ls
sleep 30 
"""

# run the job, the process ID is in out
# proc.stdin.write(job_string.encode('utf-8'))
out, err = proc.communicate()
print('out:', out)
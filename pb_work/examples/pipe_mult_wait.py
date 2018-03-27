"""=================================================================================================
run multiple commands with subprocess.Popen() and wait for all commands to finish
================================================================================================="""
import subprocess as sub
import random

jobs = []
n = 3

# start n jobs, each job sleeps a random number of up to 30 seconds, then terminates
for i in range(n):
    sec = random.randrange(30)
    command = 'sleep {}'.format(sec)
    print('job {}: {}'.format(i, command))
    job = sub.Popen(command, shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
    jobs.append(job)

# wait for all jobs to finish
print()
for i in range(n):
    jobs[i].wait()
    print('job {} finished'.format(i))
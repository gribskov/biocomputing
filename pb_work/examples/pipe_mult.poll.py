"""=================================================================================================
Run multiple jobs using subprocess.Popen()
Poll jobs until all jobs are done
================================================================================================="""
import subprocess as sub
import random
from time import sleep

jobs = []
n = 3

# start n jobs, each job sleeps a random number of seconds, then terminates
for i in range(n):
    sec = random.randrange(30)
    command = 'sleep {}'.format(sec)
    print('job {}: {}'.format(i, command))
    job = sub.Popen(command, shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
    jobs.append(job)

# poll until all jobs finish
done = 0
delay = 5  # number of seconds to wait between polls
while done < n:
    print('\nPolling')
    for i in range(n):
        if jobs[i] == 'Done':
            continue

        print('    job {} ...'.format(i), end='')
        result = jobs[i].poll()

        if result != None:
            print('finished')
            jobs[i] = 'Done'
            done += 1

        else:
            print('still running')
            # print('job {}:result={}'.format(i, result))

    sleep(delay)

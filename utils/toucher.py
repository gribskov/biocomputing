"""=================================================================================================
Walk a directory tree and update access dates in a consistent way so that time sorting is not
destroyed

8 January 2019     Michael Gribskov
================================================================================================="""
import sys
import os
import time

roots = []

if len(sys.argv) > 1:
    roots = sys.argv[1:]

else:
    # if no starting point is provided use the current working directory
    roots.append('.')

purgeday = 60
purgesec = int(0.8 * purgeday * 24 * 3600)
print( 'purge time: {} s'.format(purgesec))

for root in roots:
    for thisdir, dirlist, filelist in os.walk(root):
        print('Found directory: {}'.format(thisdir))
        now = time.time()
        oldest = now
        oldest_file = ''
        times = {}
        for file in filelist:
            # find the oldest file
            full_file = os.path.join(thisdir,file)
            times[full_file] = os.stat(full_file).st_atime
            if times[full_file] < oldest:
                oldest_file = full_file
                oldest = times[full_file]
        print('oldest file: {}\t{}\t{}\t{}'.format(oldest_file, times[oldest_file], now, int(now-times[oldest_file])))

        for file in filelist:
            print('\t{}'.format(file))
            file = os.path.join(thisdir, file)
            print('\t{}'.format(file))

            stat = os.stat(file)
            print(stat)
            print('size: {}'.format(time.ctime(stat.st_size)))
            print('accessed: {}'.format(time.ctime(stat.st_atime)))
            print('modified: {}'.format(time.ctime(stat.st_mtime)))
            print('changed: {}'.format(time.ctime(stat.st_ctime)))

            # os.utime(file, (stat.st_atime+3600, stat.st_mtime))

            stat = os.stat(file)
            print(stat)
            print('size: {}'.format(time.ctime(stat.st_size)))
            print('accessed: {}'.format(time.ctime(stat.st_atime)))
            print('modified: {}'.format(time.ctime(stat.st_mtime)))
            print('changed: {}'.format(time.ctime(stat.st_ctime)))

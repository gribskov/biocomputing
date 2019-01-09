"""=================================================================================================
Walk a directory tree and update access dates in a consistent way so that time sorting is not
destroyed

8 January 2019     Michael Gribskov
================================================================================================="""
import sys
import os
import time


def find_target_directories():
    """---------------------------------------------------------------------------------------------
    # TODO function to identify directories that need updating
    ---------------------------------------------------------------------------------------------"""
    pass
    return


roots = []

if len(sys.argv) > 1:
    roots = sys.argv[1:]

else:
    # if no starting point is provided use the current working directory
    roots.append('.')

purgeday = 60  # purge limit in days, set by RCAC
purgesec = purgeday * 24 * 3600  # purge limit in seconds
targetsec = int(
    0.5 * purgeday * 24 * 3600)  # target after access time adjustment, 50% of purge limit
print('purge limit: {} day {} sec'.format(purgeday, purgesec))
print('purge target: {} s'.format(purgesec))

for root in roots:
    for thisdir, dirlist, filelist in os.walk(root):
        print('\ndirectory: {}'.format(thisdir))
        now = time.time()
        oldest = now
        newest = 0
        oldest_file = ''
        times = {}
        for file in filelist:
            # find the oldest file
            full_file = os.path.join(thisdir, file)
            times[full_file] = os.stat(full_file).st_atime
            if times[full_file] < oldest:
                oldest_file = full_file
                oldest = times[full_file]
            if times[full_file] > newest:
                newest = times[full_file]
        # print('oldest file: {}\t{}\t{}\t{}'.format(oldest_file, times[oldest_file], now, int(now-times[oldest_file])))

        scale = (newest - oldest) / targetsec
        offset = (1 - scale) * newest

        for file in filelist:
            # print('\t{}'.format(file))
            full_file = os.path.join(thisdir, file)
            # print('\t{}'.format(file))

            stat = os.stat(full_file)
            newaccess = stat.st_atime * scale + offset

            # print(stat)
            print('accessed: {} -> {}\t{}\t{}\t\t{}\t{}'.format(
                time.ctime(stat.st_atime),
                time.ctime(newaccess),
                file,
                stat.st_size,
                time.ctime(stat.st_mtime),
                time.ctime(stat.st_ctime)))

            # TODO turn on to actually change access time
            os.utime(full_file, (stat.st_atime + 3600, stat.st_mtime))

            # stat = os.stat(full_file)
            # print(stat)
            # print('size: {}'.format(time.ctime(stat.st_size)))
            # print('accessed: {}'.format(time.ctime(stat.st_atime)))
            # print('modified: {}'.format(time.ctime(stat.st_mtime)))
            # print('changed: {}'.format(time.ctime(stat.st_ctime)))

"""=================================================================================================
Walk a directory tree and update access dates in a consistent way so that time sorting is not
destroyed. modification and creation dates are not changed.  Symbolic links are not followed, and
missing links are ignored

usage:
toucher.py <path_to_dir1> [path_to_dir2> ...]

8 January 2019     Michael Gribskov
================================================================================================="""
import sys
import os
import time


def find_target_directories(root, purgesec, old):
    """---------------------------------------------------------------------------------------------
    walk the directory hierarchy, bottom-up, finding the oldest file in the subtree under each
    directory. The list of old files is passed as an argument so multiple root directories can be
    merged into a single list
    ---------------------------------------------------------------------------------------------"""
    now = time.time()
    cutoff = now - purgesec

    for thisdir, dirlist, filelist in os.walk(root, topdown=False):
        oldest = now
        newest = 0

        need_update = False
        for file in filelist:
            try:
                # broken symbolic links fail, skip them, possibly other things?
                full_file = os.path.join(thisdir, file)
                stat = os.stat(full_file)
            except FileNotFoundError:
                continue

            if stat.st_ctime < oldest:
                oldest = stat.st_ctime
            if stat.st_ctime > newest:
                newest = stat.st_ctime
            if stat.st_atime < cutoff:
                need_update = True

        if need_update:
            old[thisdir] = (oldest, newest)

    return old


# ==================================================================================================
# main
# ==================================================================================================
roots = []

if len(sys.argv) > 1:
    roots = sys.argv[1:]

else:
    # if no starting point is provided use the current working directory
    roots.append('.')

purgeday = 46  # purge limit in days, set by RCAC, 60 days is the limit so use 60 days - 2 weeks
purgesec = purgeday * 24 * 3600  # purge limit in seconds
targetsec = int(
    0.5 * purgeday * 24 * 3600)  # target after access time adjustment, 50% of purge limit
print('purge limit: {} day {} sec'.format(purgeday, purgesec))
print('purge target: {} s'.format(targetsec))

old = {}
for root in roots:
    old = find_target_directories(root, purgesec, old)

print('{} directories containing stale files found'.format(len(old)))

for path in old:

    now = time.time()
    (oldest, newest) = old[path]
    try:
        scale = targetsec / (newest - oldest)
    except ZeroDivisionError:
        # print('\nerror: directory: {}\tnewest={}\toldest={}'.format(
        #     path, newest, oldest))
        scale = 1

    check = now - scale * (newest - oldest)
    print('\ndirectory: {}\t\t scale={:.3f}\nnewest={}\noldest={}->{}'.format(
        path, scale, time.ctime(newest), time.ctime(oldest), time.ctime(check)))

    for file in os.listdir(path):

        try:
            # failures such as broken symbolic links
            full_file = os.path.join(path, file)
            stat = os.stat(full_file)
        except FileNotFoundError:
            continue

        if os.path.isdir(full_file):
            # skip directories
            continue

        newaccess = now - scale * (newest - stat.st_ctime)

        # print(stat)
        print('accessed: {} -> {}\t{}\t{}\t\t{}\t{}'.format(
            time.ctime(stat.st_atime),
            time.ctime(newaccess),
            file,
            stat.st_size,
            time.ctime(stat.st_mtime),
            time.ctime(stat.st_ctime)))

        # change access time, leave modification time (mtime) and creation time (ctime) the same
        try:
            os.utime(full_file, (newaccess, stat.st_mtime))
        except OSError:
            sys.stderr.write('OSError:\n')

    # end of loop over files in directory
# end of loop over directories

exit(0)

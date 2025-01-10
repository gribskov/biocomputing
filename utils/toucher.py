"""=================================================================================================
Walk a directory tree and update access dates in a consistent way so that time sorting is not
destroyed. modification and creation dates are not changed.  Symbolic links are not followed, and
missing links are ignored

usage:
toucher.py <path_to_dir1> [path_to_dir2> ...]

8 January 2019     Michael Gribskov
================================================================================================="""
import argparse
import textwrap as _textwrap
import os
import sys
import time
from collections import defaultdict


# class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
class CustomFormatter(argparse.MetavarTypeHelpFormatter):
    """=============================================================================================
    Custom formatter for command line argument help
    ============================================================================================="""

    def _split_lines(self, text, width=60):
        """-----------------------------------------------------------------------------------------
        Gracefully split lines in command line help. lines are split at 60 characters by default

        :param text: str, text to split
        :param width: int, width to split at
        :return:
        -----------------------------------------------------------------------------------------"""
        text = self._whitespace_matcher.sub(' ', text).strip()
        self._max_help_position = 100
        return _textwrap.wrap(text, width)


# end of class CustomFormatter


def arguments_get():
    """---------------------------------------------------------------------------------------------
    Set up command line arguments, read from command line, and store in argparse.ArgumentParser
    object, cl. Requires customFormatter()

    :return: argparse.ArgumentParser object
    ---------------------------------------------------------------------------------------------"""

    cl = argparse.ArgumentParser(description='Update access date for old files',
                                 formatter_class=CustomFormatter,
                                 allow_abbrev=True)
    cl.add_argument('-t', '--target', type=str, default="./",
                    help='target directories (default: %(default)s)')
    cl.add_argument('-y', '--youngest_date', type=int, default=3,
                    help='Minimum age in days (default: %(default)s)')
    cl.add_argument('-o', '--oldest_date', type=int, default=21,
                    help='Maximum age in days (default: %(default)s)')
    cl.add_argument('-n', '--noupdate', action='store_true',
                    help='Report only, do not update (default: %(default)s)')

    return cl.parse_args()  # parse_args  reads the command line


def target_split(comma_str):
    """---------------------------------------------------------------------------------------------
    The target string may be a comma separated list, split on commas and return as a list

    :param comma_str: str   one or more target directories
    :return: list           with each directory as a target
    ---------------------------------------------------------------------------------------------"""
    target_list = comma_str.split(',')
    return target_list


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

args = arguments_get()
target = target_split(args.target)
ysec = args.youngest_date * 24 * 3600
osec = args.oldest_date * 24 * 3600

print(f'purge target: {args.target}')
print(f'purge range: {osec} s - {ysec} s')

# use target as a stack, visited records directories that have been seen so they don't get processed
# twice
visited = defaultdict(int)

# print('{} directories containing stale files found'.format(len(old)))

while target:
    current = target.pop()
    visited[current] += 1
    files = os.scandir(current)
    now = int(time.time())
    cutoff = now - osec
    oldest = now
    updatable = {}

    for f in files:


        if f.is_dir():
            visited[f] += 1
            target.append(f)

        # elif f.is_file:
        #     print(f'{current}\t{f.name}')

        else:
            # the other possibility is a symbolic link, do nothing
            continue

        info = os.stat(f)
        atime = int(info.st_atime)
        oldest = min(atime, oldest)
        if atime < cutoff:
            updatable[f] = atime

    scale = cutoff / oldest
    print(f'\ndirectory:{current}\tscale:{scale}\tcutoff:{cutoff}\ttoldest:{oldest}\tnow:{now}')
    for f in updatable:
        new = int(updatable[f] * scale)
        # if new > now:
        #     continue
        delta = (updatable[f] - new) / 3600 / 24
        print(f'{f}\told:{updatable[f]}\tnew:{new}\tdays:{delta}')



# for path in old:
#
#     now = time.time()
#     (oldest, newest) = old[path]
#     try:
#         scale = targetsec / (newest - oldest)
#     except ZeroDivisionError:
#         # print('\nerror: directory: {}\tnewest={}\toldest={}'.format(
#         #     path, newest, oldest))
#         scale = 1
#
#     check = now - scale * (newest - oldest)
#     print('\ndirectory: {}\t\t scale={:.3f}\nnewest={}\noldest={}->{}'.format(
#         path, scale, time.ctime(newest), time.ctime(oldest), time.ctime(check)))
#
#     for file in os.listdir(path):
#
#         try:
#             # failures such as broken symbolic links
#             full_file = os.path.join(path, file)
#             stat = os.stat(full_file)
#         except FileNotFoundError:
#             continue
#
#         if os.path.isdir(full_file):
#             # skip directories
#             continue
#
#         newaccess = now - scale * (newest - stat.st_ctime)
#
#         # print(stat)
#         print('accessed: {} -> {}\t{}\t{}\t\t{}\t{}'.format(
#             time.ctime(stat.st_atime),
#             time.ctime(newaccess),
#             file,
#             stat.st_size,
#             time.ctime(stat.st_mtime),
#             time.ctime(stat.st_ctime)))
#
#         # change access time, leave modification time (mtime) and creation time (ctime) the same
#         try:
#             os.utime(full_file, (newaccess, stat.st_mtime))
#         except OSError:
#             sys.stderr.write('OSError:\n')

    # end of loop over files in directory
# end of loop over directories

exit(0)

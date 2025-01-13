"""=================================================================================================
Starting from a comma delimited list of directory, trace the directory tree and update the access
time (st_atime) so that all times are newer than the oldest_date specified. All file and directories
are update so that the resulting files still sort the same way with respect to a_time
modification (st_mtime) and creation dates (st_ctime) are not changed.  Symbolic links are not
followed

usage: toucher.py [-h] [-t str] [-o int]

Update access date for old files

options:
  -h, --help            show this help message and exit
  -t str, --target str         target directories (default: ./)
  -o int, --oldest_date int    Maximum age in days (default: 21)

TODO implement noupdate option
TODO change directory list input into required argument

8 January 2019     Michael Gribskov
================================================================================================="""
import argparse
import textwrap as _textwrap
import os
import sys
import time
from datetime import date
from collections import defaultdict


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


# ==================================================================================================
# main
# ==================================================================================================
args = arguments_get()
target = target_split(args.target)
ysec = args.youngest_date * 24 * 3600
osec = args.oldest_date * 24 * 3600

print(f'purge target: {args.target}')
print(f'purge range: {osec} s - {ysec} s')

# use target as a stack, visited records directories that have been seen so they don't get processed
# twice. this should not happen since liknk should not be followed
visited = defaultdict(int)

while target:
    current = target.pop()
    visited[current] += 1
    files = os.scandir(current)
    now = int(time.time())
    cutoff = now - osec
    oldest = now
    updatable = {}

    # the initial directories are not DirEntry objects so current.path will fail so we need to
    # create a name for the directory currently being processed
    if isinstance(current, os.DirEntry):
        curname = current.path
    else:
        curname = current

    # processing files in curname
    count = 0
    skipped = 0
    for f in files:
        count += 1
        if f.is_dir():
            # only directories get pushed on stack, and only if they have not been visited
            visited[f] += 1
            if visited[f] <= 1:
                target.append(f)

        elif f.is_file(follow_symlinks=False):
            pass

        else:
            # the other possibility is a symbolic link, do nothing
            skipped += 1
            continue

        info = os.stat(f)
        atime = int(info.st_atime)
        oldest = min(atime, oldest)
        updatable[f] = atime

    print(
        f'\n{curname}: {len(updatable)}/{skipped} updatable/skipped. oldest file: {date.fromtimestamp(oldest)} stack={len(target)}')
    try:
        slope = (now - cutoff) / (now - oldest)
    except ZeroDivisionError:
        slope = 1

    intercept = cutoff - slope * oldest
    print(f'slope:{slope}\tintercept:{intercept}')
    n = 0

    for f in updatable:
        new = int(intercept + slope * updatable[f])
        if new > now:
            # this should never happen, but just in case
            print(f'\tWARNING:{curname}\tscale:{slope}\tcutoff:{cutoff}\toldest:{oldest}\tnow:{now}')

        # change access time, leave modification time (mtime) and creation time (ctime) the same
        try:
            n += 1
            status = os.stat(f)
            os.utime(f, (new, status.st_mtime))
            status = os.stat(f)
        except OSError:
            sys.stderr.write('OSError:\n')

exit(0)

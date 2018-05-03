import os


class Directory:
    """=============================================================================================
    Traverse a directory structure ignoring files beginning with strings listed in skipbegin and
    files files ending with strings in skipend.
    Defaults:
        self.skipbegin = ['.', '_', 'venv']
        self.skipend = ['.pyc']

    Synopsis:
        spider = Directory(dir='..\pb_work')

        while spider.next():
            print('\ndirectory: {}'.format(spider.current))
            for fname in spider.filelist:
                print('\t%s' % fname)

    31 May 2018 Michael Gribskov
    ============================================================================================="""

    def __init__(self, dir='.'):
        """-----------------------------------------------------------------------------------------
        Constructor
        -----------------------------------------------------------------------------------------"""
        self.current = None
        self.explored = {}
        self.unexplored = [dir]
        self.filelist = []
        self.file = []
        self.skipbegin = ['.', '_', 'venv']
        self.skipend = ['.pyc']

    def next(self):
        """-----------------------------------------------------------------------------------------
        Process the next directory
        :return: boolean, True if a directory is read
        -----------------------------------------------------------------------------------------"""
        if len(self.unexplored) == 0:
            return False

        self.current = self.unexplored.pop()
        self.explored[self.current] = 1
        self.filelist = []

        skip = False
        for file in os.scandir(self.current):

            # skip filenames that begin with strings in skipbegin
            skip = False
            for start in self.skipbegin:
                if file.name.startswith(start):
                    skip = True
                    break
            if skip:
                continue

            # skip filenames that end in strings in skipend
            skip = False
            for end in self.skipend:
                if file.name.endswith(end):
                    skip = True
                    break
            if skip:
                continue

            if file.is_file():
                # add files to filelist
                self.filelist.append(os.path.join(self.current, file.name))

            elif file.is_dir():
                # add directories to unexplored
                fullname = os.path.join(self.current, file.name)
                if fullname not in self.explored:
                    # prevents looping
                    self.unexplored.append(fullname)

        return True

    # end of next


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    spider = Directory(dir='..\pb_work')

    while spider.next():
        print('\ndirectory: {}'.format(spider.current))
        for fname in spider.filelist:
            print('\t%s' % fname)

    exit(0)

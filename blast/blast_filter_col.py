####################################################################################################
# check the specified column and see if value is greater than a threshold
#
####################################################################################################
import sys


def keep(fp, line):
    """---------------------------------------------------------------------------------------------

    :param line:
    :return:
    ---------------------------------------------------------------------------------------------"""
    fp.write(line)

    return


def skip(line):
    pass
    return


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    infile = sys.argv[1]
    sys.stderr.write('Blast search: {}\n'.format(infile))
    try:
        blast = open(infile, 'r')
    except OSError:
        sys.stderr.write('Cannot open blast file ({})'.format(infile))
        exit(1)

    target = 9
    threshold = 75

    # switch these to keep when column value is less
    gtfxn = keep
    ltfxn = skip

    for line in blast:

        if line.startswith('#'):
            sys.stdout.write(line)
            continue

        field = line.split("\t")
        if int(field[target]) >= threshold:
            gtfxn(sys.stdout, line)
        else:
            ltfxn(line)

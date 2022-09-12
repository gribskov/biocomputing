"""=================================================================================================
Meraculous requires the sequences to have illumina type names

/^@\S+\:\d+\:\d+\:\d+\:\d+\s+[12]\:[YN]\:\d+\:[ACTGN0]*$/
Example: @HISEQ03:379:C2WP8ACXX:7:1101:1465:2056 2:N:0:ACTTGA

Two naming conventions in the current data
1)
@SRR5992811.16 HISEQ:393:C7HT6ANXX:7:1101:2484:1913 length=125

this is easy to convert by removing the SRR tag and inventing the read1/read2 tag, e.g.
1:N:0:ACTTG
2:N:0:ACTTG

2)
@SRR5992812.7763784  7763784 length=125
@SRR5992812:run:flowcell:lane:tile:776:3784 [12]:N:0:barcode


Michael Gribskov     16 September 2020
================================================================================================="""
import sys

# import glob

# --------------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    filename = sys.argv[1]
    read = sys.argv[2]

    # defaults
    run = 111
    lane = 1
    flowcell = 'A1BC2AAXX'
    lane = 1
    tile = 1101
    barcode = 'AAAAAA'

    # for filename in glob.glob(target):
    #
    fq = None
    try:
        fq = open(filename, 'r')
    except (OSError, IOError):
        print('SRR2Illumina::open - file open error ({})'.format(filename))
        exit(1)

    nline = 0
    for line in fq:
        nline += 1
        if nline % 4 != 1:
            sys.stdout.write(line)
        else:
            field = line.rstrip().split(' ')
            if field[1].startswith('HISEQ'):
                sys.stdout.write('@{} {}:N:0:{}\n'.format(field[1], read, barcode))
            else:
                dotpos = field[0].index('.')
                readno = int(field[0][dotpos + 1:])
                y = readno % 10000
                x = readno // 10000
                field[0] = field[0][:dotpos]

                sys.stdout.write('{}:{}:{}:{}:{}:{}:{} {}:N:0:{}\n'.
                                 format(field[0], run, flowcell, lane, tile, x, y, read, barcode))

    exit(0)

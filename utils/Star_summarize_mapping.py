"""=================================================================================================
Summarize the mapping statistics in a set of STAR *.final.out files
Each file becomes a row, and each statistic a column

File example:
                                 Started job on |       Nov 29 06:34:36
                             Started mapping on |       Nov 29 06:34:37
                                    Finished on |       Nov 29 06:45:35
       Mapping speed, Million of reads per hour |       113.63

                          Number of input reads |       20769462
                      Average input read length |       300
                                    UNIQUE READS:
                   Uniquely mapped reads number |       7904955
                        Uniquely mapped reads % |       38.06%


Michael Gribskov     01 December 2022
================================================================================================="""
import sys
import os
import glob

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    pathtarget = './'
    filetarget = '*.final.out'
    target = pathtarget + '/' + filetarget

    data = {}
    for filename in glob.glob(target):
        # get the column name from the file name.  expect filenames suc as
        # C93T6R1.Log.final.out or C93T6R1.unsortedLog.final.out
        # remove path and truncate at the first period
        basename = os.path.split(filename)[1]
        period_pos = basename.find('.')
        row = basename[:period_pos]
        sys.stdout.write(f'\nfile: {filename}\trow: {row}\n')

        mapinfo = open(filename, 'r')
        data[row] = {}
        for line in mapinfo:
            if line.find('|') > -1:
                (tag, value) = line.split('|')
                print(f'{row}\t{tag}\t{value}')
                data[row][tag.strip()] = value.strip()

        mapinfo.close()

    # write out data
    summary = open('mapstats.out', 'w')
    summary.write('sample')
    for tag in data[row]:
        summary.write(f'\t{tag}')
    summary.write('\n')

    for row in data:
        summary.write(row)
        for tag in data[row]:
            summary.write('\t')
            summary.write(data[row][tag])
        summary.write('\n')


    exit(0)

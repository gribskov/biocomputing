"""=================================================================================================
Purdue genome sequencing facility data typically comprise
1) a number of directories (typically Unaligned, Unaligned_filterer, and Unaligned_filtered_unpaired
2) a number of segments per sample, 001, 002, 003, etc
3) possibly multiple lanes, L001, L002, etc.
4) possibly multiple runs, run311, run312, etc

A typical Unaligned directory might look like
-rw-r-x--- 1 mgribsko student 1183259854 Apr 30 21:30 007991_Wai-01_ATCACG_run311_L003_R1_001.fastq.gz
-rw-r-x--- 1 mgribsko student       2931 Apr 30 21:30 007991_Wai-01_ATCACG_run311_L003_R1_001.fastq.gz.stats
-rw-r-x--- 1 mgribsko student 1160272231 Apr 30 21:30 007991_Wai-01_ATCACG_run311_L003_R2_001.fastq.gz
-rw-r-x--- 1 mgribsko student       2931 Apr 30 21:30 007991_Wai-01_ATCACG_run311_L003_R2_001.fastq.gz.stats
-rw-r-x--- 1 mgribsko student  646159178 Apr 30 21:30 007991_Wai-01_ATCACG_run312_L001_R1_001.fastq.gz
-rw-r-x--- 1 mgribsko student       2903 Apr 30 21:30 007991_Wai-01_ATCACG_run312_L001_R1_001.fastq.gz.stats
-rw-r-x--- 1 mgribsko student  634029564 Apr 30 21:30 007991_Wai-01_ATCACG_run312_L001_R2_001.fastq.gz
-rw-r-x--- 1 mgribsko student       2903 Apr 30 21:30 007991_Wai-01_ATCACG_run312_L001_R2_001.fastq.gz.stats
-rw-r-x--- 1 mgribsko student       1663 Apr 30 21:30 Mycoplasma.stats.html
-rw-r-x--- 1 mgribsko student       1573 Apr 30 21:30 Mycoplasma.stats.html.all.html
-rw-r-x--- 1 mgribsko student       1255 Apr 30 21:30 Mycoplasma.stats.nohead
-rw-r-x--- 1 mgribsko student      50266 Apr 30 21:30 Mycoplasma.stats.query.summary
-rw-r-x--- 1 mgribsko student        955 Apr 30 21:30 Mycoplasma.stats.txt
drwxr-x--- 6 mgribsko student       4096 Apr 30 21:30 QC
-rw-r-x--- 1 mgribsko student       3760 Apr 30 21:30 index.html

in this case, we want to merge the two runs (same for read 2)
    007991_Wai-01_ATCACG_run311_L003_R1_001.fastq.gz
    007991_Wai-01_ATCACG_run312_L001_R1_001.fastq.gz

so we want to match on the sample: 007991_Wai-01_
and the read: R1/R2
and ignore the other fields
the files should be uncompressed and concatenated, then recompressed

31 April 2018   Michael Gribskov

================================================================================================="""
import os
import re
import subprocess
from directory import Directory


# ==================================================================================================
# main
# ==================================================================================================
if __name__ == '__main__':

    base = '../../lawson/raw'
    target_directory = re.compile('Unaligned$')
    fastq = re.compile('fastq|fq')

    spider = Directory(dir=base)
    spider.skipbegin.append('QC')
    spider.skipbegin.append('007991_Wai')

    spider.skipend.append('stats')
    spider.skipend.append('combo')
    spider.skipend.append('html')
    spider.skipend.append('json')

    while spider.next():
        if not target_directory.search(spider.current):
            continue

        print('\ncurrent:', spider.current)
        print('directory: {}'.format(spider.current))
        candidate = []
        for fname in spider.filelist:
            if fastq.search(fname):
                # print('\t%s' % fname)
                head, tail = os.path.split(fname)
                candidate.append(tail)
        candidate.sort()

        # uncompress files
        for fnum in range(len(candidate)):
            if candidate[fnum].endswith('.gz'):
                cmd = 'unpigz ' + spider.current + '/' + candidate[fnum]
                print('cmd:', cmd)
                status = subprocess.run(cmd, shell=True)
                candidate[fnum] = candidate[fnum].replace('.gz','')

        # match the candidates into groups
        idpart = re.compile('(?P<sample>.*)(?P<run>_run\d+)(?P<lane>_L\d+)(?P<readdir>_R\d)(?P<segment>_\d+)(?P<suffix>\..*)')
        # print('candidate', candidate)
        # m records the id parts of each candidate
        # set is the list of candidates in the same group (will be merged)
        m = []
        set = []

        # the first candidate is always set[0]
        set.append([0])
        m.append(idpart.search(candidate[0]))

        for f1 in range(1,len(candidate)):
            # each other candidate is checked against the existing set, if no match it creates a
            # new set
            m.append(idpart.search(candidate[f1]))
            found = False
            for s in range(len(set)):
                # the id parst of the first candidate of the set is the paradigm for the set, sm
                sm = m[set[s][0]]
                if sm['sample']==m[f1]['sample'] and sm['readdir']==m[f1]['readdir']:
                    set[s].append(f1)
                    found = True
                    break

            if not found:
                set.append([f1])

        # create commands for merging candidates by set using cat
        for s in range(len(set)):
            # print('\nmatched')
            if len(set[s]) > 1:
                cmd = 'cat '
                for fnum in set[s]:
                    id = m[fnum]['sample']
                    id += m[fnum]['run']
                    id += m[fnum]['lane']
                    id += m[fnum]['readdir']
                    id += m[fnum]['segment']
                    id += m[fnum]['suffix']
                    # print('    ', id)
                    cmd += spider.current + '/' + id + ' '
                lead = set[s][0]
                out = spider.current + '/' + m[lead]['sample']+m[lead]['readdir']+'_merged'+m[lead]['suffix']
                cmd += ' > ' + out
                print('cmd:', cmd)
                status = subprocess.run(cmd, shell=True)


        # recompress all fastq/fq files in directory
        cmd1 = 'pigz ' + spider.current + '/*.fastq'
        cmd2 = 'pigz ' + spider.current + '/*.fq'
        print('cmd:', cmd1)
        print('cmd:', cmd2)
        status = subprocess.run(cmd1, shell=True)
        status = subprocess.run(cmd2, shell=True)


    exit(0)


#!/bin/env python3
"""-------------------------------------------------------------------------------------------------
The NCBI naming convention for reads does not allow Trinity to identify read 1 and read 2.
this script renames the read using the trinity convention /1 and /2 for reads 1 and 2, respectively.
All files matching the file name protorype, e.g. SRR*.paired.fastq, will be converted and output
placed in, e.g. SRR*.paired.renamed.fastq.

usage:
    srr_to_trinity.py <file_name_prototype>

Michael Gribskov    7   August 2018
-------------------------------------------------------------------------------------------------"""
import sys
import os
import glob

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    sys.stderr.write('files: {}\n'.format(sys.argv[1]))

    for wild in glob.iglob(sys.argv[1]):
        path, file = os.path.split(wild)
        sys.stderr.write('processing file: {}\n'.format(file))

exit(0)

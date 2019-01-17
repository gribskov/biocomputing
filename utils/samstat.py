"""=================================================================================================
samstat files are organized one per file and by row, this script reorganizes into sequences in
columns

SN      raw total sequences:    85128836
SN      filtered sequences:     0
SN      sequences:      85128836
SN      is sorted:      1
SN      1st fragments:  42564418
SN      last fragments: 42564418
SN      reads mapped:   84909236
SN      reads mapped and paired:        84768496        # paired-end technology bit set + both mates mapped
SN      reads unmapped: 219600
SN      reads properly paired:  83828528        # proper-pair bit set
SN      reads paired:   85128836        # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      215213  # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 0
SN      total length:   12854454236     # ignores clipping
SN      bases mapped:   12821294636     # ignores clipping
SN      bases mapped (cigar):   12093112251     # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     78411927        # from NM fields
SN      error rate:     6.484015e-03    # mismatches / bases mapped (cigar)
SN      average length: 151
SN      maximum length: 151
SN      average quality:        37.7
SN      insert size average:    560.6
SN      insert size standard deviation: 209.4
SN      inward oriented pairs:  41895271
SN      outward oriented pairs: 79993
SN      pairs with other orientation:   68834
SN      pairs on different chromosomes: 284295

Michael Gribskov     16 January 2019
================================================================================================="""
import os

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # find all matching files in the current directory
    suffix = '.stat'

    for file in os.scandir():
        if file.name.endswith(suffix):
            print('yes', file)
        else:
            print('no', file)


exit(0)
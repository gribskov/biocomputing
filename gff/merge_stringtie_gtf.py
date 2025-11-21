"""=================================================================================================
Merge overlapping transcripts in GTF file to produce transcripts suitable for use with Salmon
(or Kallisto). Works on a single input file such as the one procuced by stringtie merge

1. Only records with the desired feature(s), usually gene or transcript
2. merged features begin at the earliest beginning and end at the latest end point
================================================================================================="""
from gff import Gff
from ranges.linear_range import Range

# ##################################################################################################
# main program
# ##################################################################################################
if __name__ == '__main__':
    gtffile = 'data/stringtie_merged_jm.gtf'
    features = ['transcript']

    print(f'merge_stringtie.py')
    print(f'\tGTF file: {gtffile}')
    print(f'\tSelected features: {features}')

    # read in gtf with selected features
    gtf = Gff(file=gtffile)
    merged = Range()
    feature_n = 0
    for f in gtf.get_by_feature('transcript'):
        feature_n += 1
        merged.features.append(f)



    # print(f'\t{features} features read: {len(gtf.data)}')

    # overlap

    # created merged names and write out as gtf

    exit(0)

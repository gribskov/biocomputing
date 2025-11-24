"""=================================================================================================
Merge overlapping transcripts in GTF file to produce transcripts suitable for use with Salmon
(or Kallisto). Works on a single input file such as the one produced by stringtie merge

1. Only records with the desired feature(s), usually gene or transcript
2. merged features begin at the earliest beginning and end at the latest end point
================================================================================================="""
from gff2 import GxfSet
from ranges.linear_range import Range


class Lrange:
    """---------------------------------------------------------------------------------------------

    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        self.name = ''
        self.seqid = ''
        self.start = None
        self.end = None
        self.strand = None
        self.members = []


def overlap(flist, space=1000):
    """---------------------------------------------------------------------------------------------
    find overlap between features
    
    :param features: GxfSet     collection of GxfRecord features
    :return: list               merged features, collection of Lrange
    ---------------------------------------------------------------------------------------------"""
    merged = []
    current = Lrange()

    range_n = 0
    feature_n = 0
    for f in sorted(flist.features, key=lambda x: (x.seqid, x.strand, x.start)):
        # if f.seqid != 'LG01':
        #     continue

        feature_n += 1
        if (f.seqid != current.seqid or
                f.strand != current.strand or
                f.start - current.end > space):

            # start new range
            if current.name:
                merged.append(current)
            current = Lrange()
            current.name = f'range_{range_n}'
            range_n += 1
            current.seqid = f.seqid
            current.start = f.start
            current.end = f.end
            current.strand = f.strand
            current.members.append(f)
            continue

        # merge
        current.end = max(current.end, f.end)
        current.members.append(f)

        print(f)

    return merged


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
    gtf = GxfSet(file=gtffile)
    feature_n = gtf.feature_get(features)
    print(f'{feature_n} features read from {gtffile}')

    # overlap
    merged = overlap(gtf)
    print(f'merged: {len(merged)}')

    # create merged names and write out as gtf

    output = 'merged.out'
    mout = open(output, 'w')
    for f in sorted(merged, key=lambda x:(x.seqid, x.start)):
        attr_out = {'gene_id':f.members[0].attribute['gene_id'], 'transcript_id':f.members[0].attribute['transcript_id']}
        mlist = []
        for ff in f.members:
            mlist.append(ff.attribute['transcript_id'])
        attr_out['attribute'] = ','.join(mlist)
        attr_str = '; '.join(attr_out)
        print()

        mout.write(f'{f.seqid}\tStringTie\ttranscript\t{f.start}\t{f.end}\t.\t{f.strand}\t.\t{attr_str}\n')

    exit(0)

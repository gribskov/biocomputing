"""#################################################################################################
fastq_split
break a fastq file into n pieces
#################################################################################################"""
import sys
import os
import glob
import re

"""#################################################################################################
main
#################################################################################################"""
if __name__ == '__main__':
    segments = int(sys.argv[1])
    infile = sys.argv[2]

    sys.stderr.write('input:{}\n'.format(infile))
    sys.stderr.write('segments:{}\n'.format(segments))

    for fq in glob.iglob(infile):
        base = os.path.basename(fq)
        base = re.sub('\.fq$|\.fastq$', '', base)
        sys.stderr.write('base:{}\n'.format(base))
        fqin = open(fq, 'r')

        fqout = []
        for i in range(segments):
            outfile = '{}.{:03d}.fq'.format(base, i)
            sys.stderr.write('    out[{}]: {}\n'.format(i, outfile))
            fqout.append(open(outfile,'w'))

        nline = 0
        entry = ''
        seq = 0
        for line in fqin:
            nline += 1
            entry += line
            if not nline % 4:
                fqout[seq].write('{}'.format(entry))
                entry  = ''
                seq = ( seq + 1)  % segments

            # if nline == 400:
            #     break
        for f in fqout:
            f.close()


exit(0)

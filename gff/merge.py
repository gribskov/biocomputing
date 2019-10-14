from gff import Gff

if __name__ == '__main__':

    gff = Gff(file="stringtie.gff")
    transcript_n = gff.read_feature('transcript')
    print('{} features read'.format(transcript_n))

    bundle_n = 0
    bundle = []
    b_start = 0
    b_stop = 0
    for t in sorted(gff.data, key=lambda k: int(k['begin'])):
        id = t['transcript_id']
        begin = int(t['begin'])
        end = int(t['end'])
        strand = t['strand']

        if begin > b_stop:
            # new bundle
            bundle.append([t])

            print('bundle {}'.format(bundle_n))
            bundle_n += 1
            b_start = begin
            b_stop = end

        else:
            # continue current bundle
            b_stop = end
            bundle[-1].append([t])

        print('\t{}\t{}\t{}\t{}'.format(id, begin, end, strand))

    exit(0)

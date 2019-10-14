from gff import Gff

class Bundle:

    def __init__(self):
        """

        """
        self.begin = 0
        self.end = 0
        self.sequence = ''
        self.transcript = []

def make_bundle(gff):
    """---------------------------------------------------------------------------------------------

    :param gff: Gff object
    :return: list of bundles
    ---------------------------------------------------------------------------------------------"""
    bundle_list = []
    b_stop = 0
    for t in sorted(gff.data, key=lambda k: (k['sequence'], int(k['begin']))):
        id = t['transcript_id']
        begin = int(t['begin'])
        end = int(t['end'])
        sequence = t['sequence']
        strand = t['strand']

        if begin > b_stop:
            # new bundle
            b = Bundle()
            b.begin = begin
            b.end = end
            b.sequence = sequence
            b.transcript.append(t)
            bundle_list.append(b)

        else:
            # continue current bundle
            b.stop = end
            b.transcript.append(t)

        b_stop = end

    return bundle_list

if __name__ == '__main__':

    gff = Gff(file="stringtie.gff")
    transcript_n = gff.read_feature('transcript')
    print('{} features read'.format(transcript_n))
    gff.replace_columns_re(['sequence'], 'lcl\|', r'')

    bundle = make_bundle(gff)
    bundle_n = 0
    print('{} bundles found'.format(len(bundle)))
    for b in bundle:
        print('\nbundle {}'.format(bundle_n))
        bundle_n += 1
        for transcript in b.transcript:
            id = transcript['transcript_id']
            begin = int(transcript['begin'])
            end = int(transcript['end'])
            sequence = transcript['sequence']
            strand = transcript['strand']
            print('\t{}\t{}\t{}\t{}\t{}'.format(id, transcript['sequence'], transcript['begin'], transcript['end'], transcript['strand']))

    exit(0)

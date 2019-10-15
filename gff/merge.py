from gff import Gff

class Bundle:

    def __init__(self):
        """

        """
        self.begin = 0
        self.end = 0
        self.sequence = ''
        self.transcript = []

def make_bundle(gff,id_column='transcript_id'):
    """---------------------------------------------------------------------------------------------

    :param gff: Gff object
    :param id_column: name of id column in gff hash
    :return: list of bundles
    ---------------------------------------------------------------------------------------------"""
    bundle_list = []
    b_stop = 0
    sequence_old = ''
    for t in sorted(gff.data, key=lambda k: (k['sequence'], int(k['begin']))):
        id = t[id_column]
        begin = int(t['begin'])
        end = int(t['end'])
        sequence = t['sequence']
        strand = t['strand']

        if begin > b_stop or sequence != sequence_old:
            # new bundle
            b = Bundle()
            b.begin = begin
            b.end = max(b.end,end)
            b.sequence = sequence
            b.transcript.append(t)
            bundle_list.append(b)
            sequence_old = sequence

        else:
            # continue current bundle
            b.stop = end
            b.transcript.append(t)

        b_stop = end

    return bundle_list

if __name__ == '__main__':

    ugff = Gff(file="stringtie.gff")
    transcript_n = ugff.read_feature('transcript')
    print('{} unstranded features read'.format(transcript_n))
    ugff.replace_columns_re(['sequence'], 'lcl\|', r'')
    ubundle = make_bundle(ugff)

    sgff = Gff(file="genome.gff")
    sgff.attr_sep = '='
    transcript_n = sgff.read_feature('mRNA')
    print('{} stranded features read'.format(transcript_n))

    # query = r'(maker|augustus|masked|processed|gene)(-|_)'
    query = r'(maker|augustus|masked|processed|gene|trnascan)(-|_)'

    sgff.replace_columns_re(['Parent', 'ID'], query, r'')
    sgff.replace_columns_re(['Parent', 'ID'], r'-tRNA-', r'.tRNA')
    sgff.rename_key('ID', 'transcript_id')
    sbundle = make_bundle(sgff)
    # sbundle = make_bundle(sgff, id_column='ID')


    set = []
    set_end = 0
    set_seq = ''
    for b in sorted(ubundle+sbundle, key=lambda k: (k.sequence, k.begin)):

        if b.begin > set_end or b.sequence != set_seq:
            # new set
            set_seq = b.sequence
            set_end = b.end
            set.append( {'begin'   : b.begin,
                        'end'     : set_end,
                        'sequence': set_seq,
                        'member'  : [b]} )

        else:
        # continue current set
            set[-1]['end'] = max(set[-1]['end'], b.end)
            set[-1]['member'].append(b)
            set_end = set[-1]['end']

    s_n = 0
    for s in set:
        print('set {}'.format(s_n))
        s_n += 1
        for b in s['member']:
            print('\t{}\t{}\t{}'.format(b.sequence, b.begin, b.end))
            for transcript in b.transcript:
                    id = transcript['transcript_id']
                    begin = int(transcript['begin'])
                    end = int(transcript['end'])
                    sequence = transcript['sequence']
                    strand = transcript['strand']
                    print('\t \t{}\t{}\t{}\t{}\t{}'.format(id, transcript['sequence'], transcript['begin'], transcript['end'], transcript['strand']))

    # bundle_n = 0
    # print('{} bundles found'.format(len(bundle)))
    # for b in bundle:
    #     print('\nbundle {}'.format(bundle_n))
    #     bundle_n += 1
    #     for transcript in b.transcript:
    #         id = transcript['transcript_id']
    #         begin = int(transcript['begin'])
    #         end = int(transcript['end'])
    #         sequence = transcript['sequence']
    #         strand = transcript['strand']
    #         print('\t{}\t{}\t{}\t{}\t{}'.format(id, transcript['sequence'], transcript['begin'], transcript['end'], transcript['strand']))

    exit(0)

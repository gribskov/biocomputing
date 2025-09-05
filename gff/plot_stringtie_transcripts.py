"""=================================================================================================
plot transcripts from GTF file. This is version is specific for stringtie merged GTF files

Michael Gribskov     05 September 2025
================================================================================================="""
from gff import Gff
import matplotlib.pylab as plt

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gtf = 'data/merged.gtf'
    target_feature = ['transcript', 'exon']
    sequence = ['LG01']
    selector = {'sequence': ['LG01'],
                'feature':  ['transcript', 'exon'],
                'end':      50000}
    groupby = 'gene_id'

    anno = Gff(file=gtf)
    anno.del_attr = True
    n_selected = anno.get_gene_id_groups_by_generic_selector(selector, groupby)
    print(f'{n_selected} regions selected')

    # fig, ax = plt.subplots()
    #
    # y = -2
    # current = ''
    # width = []
    # left = []
    # ypos = []
    # label = []
    # color = []
    #
    # for feature in anno.data:
    #     # plot desired features in groups of transcipt_id with all exons on a
    #     # single line
    #     if feature['transcript_id'] != current:
    #         if current:
    #             ax.barh(ypos, width=width, left=left, color=color)
    #         current = feature['transcript_id']
    #         width = []
    #         left = []
    #         ypos = []
    #         color = []
    #         # label = []
    #         y += 2
    #
    #     print(f'{feature['sequence']}\t{feature['feature']}\t{feature['transcript_id']}\t{feature['begin']}', end='')
    #     print(f'\t{feature['end']}\t{feature['strand']}')
    #     width.append(feature['end'] - feature['begin'] + 1)
    #     left.append(feature['begin'])
    #
    #     if feature['feature'] == 'transcript':
    #         ypos.append(y + 1)
    #         label.append('')
    #         color.append('0.75')
    #         label.append(feature['transcript_id'])
    #     else:
    #         ypos.append(y)
    #         color.append('r')
    #
    # ax.barh(ypos, width=width, left=left, color=color)
    # plt.yticks(range(len(label)), labels=label)
    #
    # ax.set(xlabel='Position', ylabel='Transcript', title='Exon and Transcript position')
    # ax.grid()
    #
    # # fig.savefig("test.png")
    # plt.show()

    exit(0)

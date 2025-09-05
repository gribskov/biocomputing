"""=================================================================================================
plot transcripts from GTF file. This is a pretty generic version that should be modifiable for a
specific purpose

Michael Gribskov     04 September 2025
================================================================================================="""
from gff import Gff
import matplotlib.pylab as plt


def split_attributes_gtf(anno):
    """---------------------------------------------------------------------------------------------
    Split attributes in GTF format. Slit string on ;, then extract tag-value pairs

    :param anno: Gff object
    :return: dict               dictionary of tags and values
    ---------------------------------------------------------------------------------------------"""
    for feature in anno.data:
        field = feature['attribute'].split(';')
        for f in field:
            tag, field = f.split(' ', maxsplit=1)
            feature[tag] = field.strip('"')

    return True


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gtf = 'data/merged.gtf'
    begin = 1
    end = 80000
    target_feature = ['transcript', 'exon']
    sequence = ['LG01']

    anno = Gff(file=gtf)
    anno.del_attr = True
    anno.get_by_sequence_range(target_feature, sequence, begin, end)
    for feature in anno.data:
        print(feature)

    fig, ax = plt.subplots()

    y = -2
    current = ''
    width = []
    left = []
    ypos = []
    label = []
    color = []

    for feature in anno.data:
        # plot desired features in groups of transcipt_id with all exons on a
        # single line
        if feature['transcript_id'] != current:
            if current:
                ax.barh(ypos, width=width, left=left, color=color)
            current = feature['transcript_id']
            width = []
            left = []
            ypos = []
            color = []
            # label = []
            y += 2

        print(f'{feature['sequence']}\t{feature['feature']}\t{feature['transcript_id']}\t{feature['begin']}', end='')
        print(f'\t{feature['end']}\t{feature['strand']}')
        width.append(feature['end'] - feature['begin'] + 1)
        left.append(feature['begin'])

        if feature['feature'] == 'transcript':
            ypos.append(y + 1)
            label.append('')
            color.append('0.75')
            label.append(feature['transcript_id'])
        else:
            ypos.append(y)
            color.append('r')

    ax.barh(ypos, width=width, left=left, color=color)
    plt.yticks(range(len(label)), labels=label)

    ax.set(xlabel='Position', ylabel='Transcript', title='Exon and Transcript position')
    ax.grid()

    # fig.savefig("test.png")
    plt.show()

    exit(0)

"""=================================================================================================
plot transcripts from GTF file

Michael Gribskov     04 September 2025
================================================================================================="""
from gff import Gff


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
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gtf = 'data/merged.gtf'
    anno = Gff(file=gtf)
    anno.del_attr = True
    anno.get_by_sequence_range(['transcript', 'exon'], ['LG01'], 1, 200000)
    for feature in anno.data:
        print(feature)

    for feature in sorted(anno.data, key=lambda x: x['begin']):
        print(f'{feature['sequence']}\t{feature['feature']}\t{feature['begin']}\t{feature['end']}')

    exit(0)

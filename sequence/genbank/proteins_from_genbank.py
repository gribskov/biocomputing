"""=================================================================================================
extract protein coding sequences from genbank nucleotide file.  Designed for mitochondrial genome

Michael Gribskov     15 October 2022
================================================================================================="""
import sys


def get_id(file, taxa_names):
    """---------------------------------------------------------------------------------------------
    get a unique ID from the description line

    :param taxa_names: list     names of taxa
    :param file: filehandle     source genbank file
    :return: logical            True if id could be read, false at EOF
    ---------------------------------------------------------------------------------------------"""
    remove = 'isolate breed voucher strain: mitochondrion complete partial genome almost mitochondrial Hybrid horse historical'
    # punc = {'.': '', ',': '', '\n':''}
    punc = ''.maketrans('\\/:_', '___ ', '.,())')

    desc = ''
    desc_found = None
    for line in file:
        if line.startswith('DEFINITION'):
            line = line.replace('DEFINITION ', '').rstrip()
            line = line.replace('DNA', '')
            # print(f'\n{line.rstrip()}')
            desc = line
            desc_found = True
            break
        desc_found = False

    filtered = []
    if desc_found:
        # desc = desc.replace('_', ' ')
        words = desc.translate(punc).split()
        for w in words:
            if w not in remove:
                filtered.append(w)

    id = '_'.join(filtered)
    id = id.replace('haplogroup_', 'hap')
    id = id.replace('haplotype_', 'hap')

    if id in taxa_names:
        base, n = id.rsplit('_', maxsplit=1)
        if n.isdigit():
            id = f'{base}_{int(n) + 1}'
        else:
            id = f'{id}_1'

    print(f'{id}\t\t{line}')
    taxa_names.append(id)

    return desc_found


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gb = open(sys.argv[1], 'r')
    taxa_names = []

    while get_id(gb, taxa_names):
        pass

    exit(0)

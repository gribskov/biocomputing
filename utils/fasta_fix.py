"""=================================================================================================
fix fasta files:
    rename duplicates to unique names

Michael Gribskov     23 April 2022
================================================================================================="""
import sys
# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    fasta_in = open(sys.argv[1], 'r')
    ulist = set()
    n = 0
    n_dup = 0
    for line in fasta_in:
        if line.startswith('>'):
            if line.find('GNL|')>-1:
                print(line)

            n += 1
            try:
                id, doc = line.split(' ', maxsplit=1)
            except:
                id = line
                doc = ''

            # print(f'id:{id}\tdoc:{doc}')
            # print(f'\tid:{id}')

            if id in ulist:
                print(f'{n}\t{id} is duplicate')
                n_dup += 1
            else:
                ulist.add(id)

            if not n%100000:
                print(f'{n}')

    print(f'{n} sequences found')
    print(f'{n_dup} duplicates found')




    exit(0)

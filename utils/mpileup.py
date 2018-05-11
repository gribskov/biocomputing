"""=================================================================================================
process mpileup into a simpler to read format

9 May 2018  Michael Gribskov
================================================================================================="""
import sys
import copy

def count(seq):
    """---------------------------------------------------------------------------------------------
    Count the bases in the sequence. convert lowercase, which means pposite strand to uppercase
    ignore special mpileup charcters such as $, ^, etc.

    :param seq:
    :return: dict, key=base, value=count
    ---------------------------------------------------------------------------------------------"""

    count = {}
    for base in seq:
        if base not in 'AaCcGgTt':
            continue

        try:
            count[base.upper()] += 1
        except KeyError:
            count[base.upper()] = 1

    return count

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    try:
        mp = open(sys.argv[1], 'r')
    except IndexError:
        print(' you must enter a filename on the command line')
    except FileNotFoundError:
        print(' can\'t find the specified file ({})'.format(sys.argv[1]))
    except Exception as e:
        print('unexpected exception: {}'.format(e))

    threshold = 0.95
    print('SNP threshold =', threshold)
    print('\n\t\tpos\tdepth\tfrac\tbases')
    fasta = ''
    snp = []
    for line in mp:
        # print(line)
        field = line.split()
        # fields are id, pos, ref_base, depth, aligned bases, aligned qual
        pos = int(field[1])
        basecount = count(field[4])
        depth = int(field[3])
        # print('{}\t{}\t{}'.format(field[1], depth, basecount))
        if depth:
            maxbase = max(basecount, key=basecount.get)
            frac = basecount[maxbase] / depth
            if frac < threshold:
                print('snp:\t{}\t{}\t{:.3f}\t{}'.format(pos, depth, frac, basecount))
                thiscount = copy.deepcopy(basecount)
                snp.append([pos,thiscount])
            fasta += maxbase
        else:
            fasta += 'N'

    # for s in snp:
    #     print(s)

    pos = 0
    linelen = 100
    while pos < len(fasta):
        print(fasta[pos:pos+linelen])
        pos += linelen

    exit(0)


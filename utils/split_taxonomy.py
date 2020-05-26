"""=================================================================================================
For a typical metagenomics taxonomy file like

D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Marine Benthic Group D and DHVEG-1;__
D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Marine Group II;D_4__uncultured archaeon
D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Marine Group II;D_4__uncultured euryarchaeote
D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Marine Group II;D_4__uncultured haloarchaeon
D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Marine Group II;__
D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;__;__
D_0__Archaea;D_1__Nanoarchaeaeota;D_2__Woesearchaeia;D_3__Candidatus Pacearchaeota archaeon CG1_02_32_21;D_4__Candidatus Pacearchaeota archaeon CG1_02_32_21
D_0__Archaea;D_1__Nanoarchaeaeota;D_2__Woesearchaeia;D_3__Candidatus Staskawiczbacteria bacterium RIFOXYA2_FULL_32_7;D_4__Candidatus Staskawiczbacteria bacterium RIFOXYA2_FULL_32_7
D_0__Archaea;D_1__Nanoarchaeaeota;D_2__Woesearchaeia;__;__

1. split on taxonomic divisions (semicolon)
2. write out two columns, tab-delimited
    a. nl, of the original taxonomic divisions
    b. nr, one level from the original taxonomy with
       spaces deleted
       names like __ reconstructed to be the last known level plus unknown


Michael Gribskov     21 May 2020
================================================================================================="""
import sys

def getlevel( level, target ):
    """---------------------------------------------------------------------------------------------
    from the list of levels extract the desired target level
    1) remove the level indicator, such as D_3__
    2) if the level is blank, indicated as '__', search levels to the left until a defined level is
       found.  Use this as the name with "unknown" appended
    3) remove all spaces

    :param level: list, split taxonomy string
    :param: target, int: level to select as a group label
    :return: string
    ---------------------------------------------------------------------------------------------"""

    label = ''
    depth = target
    while not label:
        token = level[depth-1].split('__')
        label = '_'.join(token[1:])
        depth -= 1

    if depth != target-1:
        label = 'unknown' + label

    return label.replace(' ','')

def level2taxonomy(level, nl):
    """---------------------------------------------------------------------------------------------
    Prune the taxonomy at the indicated level nl

    :param level: list, split taxonomy string
    :param nl: int, number of levels in final list
    :return: list, retained levels
    ---------------------------------------------------------------------------------------------"""
    if nl == 0:
        return level

    while len(level) > nl:
        level.pop()

    return level

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    nl = 0  # zero indicates all
    nr = 3

    # input taxonomy is argument 1
    tax_infile = sys.argv[1]
    sys.stderr.write('\tTaxonomy: {}\n'.format(tax_infile))
    try:
        taxon = open(tax_infile, 'r')
    except Exception as e:
        sys.stderr.write('\nUnable to open taxonomy ({})\n'.format(tax_infile))
        sys.stderr.write('{}\n'.format(e))
        exit(1)

    # output file is same name with suffix .meta
    tax_outfile = tax_infile + '.meta'
    sys.stderr.write('\tMetadata: {}\n'.format(tax_outfile))
    try:
        meta = open(tax_outfile, 'w')
    except Exception as e:
        sys.stderr.write('\nUnable to open taxonomy metadata output({})\n'.format(tax_outfile))
        sys.stderr.write('{}\n'.format(e))
        exit(1)

    meta.write('#OTU ID\tgroup\n')
    nline = 0
    for line in taxon:
        line = line.strip()
        if line.startswith('#'):
            continue;
        if not line:
            continue

        nline += 1
        level = line.split(';')
        #sys.stderr.write('{}{}\n'.format(nline, level))
        level = level2taxonomy(level, nl)
        group = getlevel(level, nr)

        meta.write('{}\t{}\n'.format(';'.join(level), group))

    sys.stderr.write('\n{} tax read from {}\n'.format(nline, tax_infile))

    # final clean up
    taxon.close()
    meta.close()

    exit(0)
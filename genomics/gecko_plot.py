"""=================================================================================================
plot the csv output file from gecko comparison of two genomes

example:
Total CSB: 0
========================================================
Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY
Frag,2219163,24003925,2220409,24005171,f,0,1247,4892,1235,98.08,99.04,0,22
Frag,2445521,24227503,2448339,24230321,f,0,2819,11084,2795,98.30,99.15,0,22
Frag,338024,22108613,339668,22110257,f,0,1645,6452,1629,98.05,99.03,0,22

Michael Gribskov     20 May 2024
================================================================================================="""
import sys
import matplotlib.pylab as plt
import matplotlib.collections as mc


def gecko_read(filename):
    """---------------------------------------------------------------------------------------------
    Read in the comparison result, see program header
    
    :param filename: 
    :return: 
    ---------------------------------------------------------------------------------------------"""
    try:
        fp = open(filename, 'r')
    except OSError:
        sys.stderr.write(f'gecko_read: Error opening input file ({filename}\n')
        exit(1)

    # skip first three lines
    fp.readline()
    fp.readline()
    fp.readline()

    match = []
    for line in fp:
        field = line.rstrip().split(',')
        row = {}
        row['type'] = field[0]
        row['xstart'] = int(field[1])
        row['ystart'] = int(field[2])
        row['xend'] = int(field[3])
        row['yend'] = int(field[4])
        row['strand(f/r)'] = field[5]
        row['block'] = int(field[6])
        row['length'] = int(field[7])
        row['score'] = int(field[8])
        row['ident'] = field[9]
        row['similarity'] = float(field[10])
        row['%ident'] = float(field[11])
        row['seqx'] = int(field[12])
        row['seqy'] = int(field[13])
        match.append(row)

    return match


def getlines(match):
    """---------------------------------------------------------------------------------------------
    make a collection of line segments for plotting. each segment is (xstart,xend) to (ystart,yend)
    
    :param match: list of dict      gecko data from gecko_read
    :return: list                   [[[xstart,xend],[ystart,yend]], ...]
             list                   xrange [xmin, xmax]
             list                   yrange [ymin, ymax]
    ---------------------------------------------------------------------------------------------"""
    xrange = [10000000000, 0]
    yrange = [10000000000, 0]
    seglist = []
    for m in match:
        seg = [[m['xstart'], m['xend']], [m['ystart'], m['yend']]]
        seglist.append(seg)
        xrange[0] = min(xrange[0], m['xstart'])
        xrange[1] = max(xrange[1], m['xend'])
        yrange[0] = min(yrange[0], m['ystart'])
        yrange[1] = max(yrange[1], m['yend'])

    return seglist, xrange, yrange


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    plot = True
    reorder = False

    geckofile = 'data/scaffolds.reordered-GCA_000814965.1_MBR_1.0_genomic.csv'
    match = gecko_read(geckofile)
    print(f'{len(match)} lines read from geckofile')
    lines, xrange, yrange = getlines(match)

    if reorder:

        # xseqfile = 'data/GCA_000814965.1_MBR_1.0_genomic.fna'
        xseqfile = 'data/scaffolds.fasta'
        xseq = open(xseqfile, 'r')
        xordfile = xseqfile.replace('fasta', 'reordered.fa')
        xord = open(xordfile, 'w')

        # read sequence file
        contigs = []
        contig_n = 0
        base_n = 0
        for line in xseq:
            if line.startswith('>'):
                # docline
                try:
                    id, doc = line.split(' ', maxsplit=1)
                except ValueError:
                    id = line
                    doc = ''
                contigs.append({'id': id, 'doc': doc, 'seq': ''})
                seq = contigs[-1]['seq']
                contig_n += 1
            else:
                contigs[-1]['seq'] += line
                base_n += len(line) - 1

        print(f'xcontigs read from {xseqfile}: {contig_n} contigs\t {base_n} bases')

        # sort by xstart
        o = open('match.out', 'w')
        xused = []
        yused = []
        for m in sorted(match, key=lambda m: m['ystart']):
            o.write(f'{m["seqx"]}\t{m["seqy"]}\t{m["xstart"]}\t{m["ystart"]}\t{m["length"]}\n')
            print(f'{m["seqx"]}\t{m["seqy"]}\t{m["xstart"]}\t{m["ystart"]}\t{m["length"]}')
            if m["seqx"] not in xused:
                xused.append(m["seqx"])
                xord.write(contigs[m["seqx"]]["id"])
                xord.write(contigs[m["seqx"]]["doc"])
                xord.write(contigs[m["seqx"]]["seq"])

            if m["seqy"] not in yused:
                yused.append(m["seqy"])

        print(f'x sequences used: {len(xused)}')
        print(f'y sequences used: {len(yused)}')

    if plot:
        fig, ax = plt.subplots()
        n = 0
        for m in match:
            ax.plot([m['xstart'], m['xend']], [m['ystart'], m['yend']])
            # n += 1
            # if n > 2:
            #     break

        ax.set(xlabel='Genome X', ylabel='Genome Y',
               title='gecko matches')
        ax.grid()

        fig.savefig("test.png")
        plt.show()

    exit(0)

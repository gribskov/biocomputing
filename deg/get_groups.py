"""=================================================================================================
groups of differentially expressed genes from R are used to select genes from the merged file
the gene file, written by write.txble(),  should look like
"x"
"1"	"TRINITY_DN21491_c0_g2_i1"
"2"	"TRINITY_DN5632_c0_g3_i5"
"3"	"TRINITY_DN4706_c0_g2_i7"
"4"	"TRINITY_DN4724_c0_g3_i1"
"5"	"TRINITY_DN4776_c0_g2_i3"
"6"	"TRINITY_DN4729_c0_g1_i2"
"7"	"TRINITY_DN4765_c0_g3_i1"
"8"	"TRINITY_DN4793_c4_g1_i2"
"9"	"TRINITY_DN4736_c0_g2_i8"
"10"	"TRINITY_DN4738_c0_g1_i3"
"11"	"TRINITY_DN4740_c0_g1_i3"
"12"	"TRINITY_DN19789_c0_g4_i2"
"13"	"TRINITY_DN44_c0_g1_i15"
"14"	"TRINITY_DN32_c0_g1_i7"
"15"	"TRINITY_DN10_c0_g1_i3"
"16"	"TRINITY_DN63_c0_g1_i4"

wher each group begins with "x"
================================================================================================="""
import sys

groupfile = sys.argv[1]
infofile = sys.argv[2]

group = None
try:
    group = open(groupfile, 'r')
except OSError:
    sys.stderr.write('Unable to open groups file ({})'.format(groupfile))
    exit(1)

try:
    info = open(infofile, 'r' )
except OSError:
    sys.stderr.write('Unable to open information file ({})'.format(infofile))

group_ids = []
ngroup = 0
for line in group:
    line = line.rstrip()
    if line == '"x"':
        # begin reading group
        group_ids.append([])
        ngroup += 1
    else:
        line = line.replace('"', '')
        tab = line.find('\t')
        group_ids[ngroup - 1].append(line[tab+1:])

for i in range(len(group_ids)):
    sys.stdout.write('group {} - {} transcripts\n'.format(i+1, len(group_ids[i])))

group.close()

# open a file for each group
gf = []
for i in range(len(group_ids)):
    filename = 'group{}'.format(i)
    try:
        gf.append(open(filename, 'w'))
    except OSError:
        sys.stderr.write('Unable to open group output file ({})'.format(filename))

bigline = ''
for line in info:

    if line.startswith('TRINITY'):
        if bigline:

            id, data = bigline.split('\t', maxsplit=1)
            # print(id, ':', data)
            for g in range(len(group_ids)):
                if id in group_ids[g]:
                    print('id {}\t group {}'.format(id, g))
                    gf[g].write('{}\t{}'.format(id,data))
                    break

        bigline = line
    else:
        bigline += line

if bigline:

    id, data = bigline.split('\t', maxsplit=1)
    # print(id, ':', data)
    for g in range(len(group_ids)):
        if id in group_ids[g]:
            print('id {}\t group {}'.format(id, g))
            gf[g].write('{}\t{}'.format(id,data))
            break

####################################################################################################
# compare transcriptome to genome.  for each transcript find begin/end of good matching region in
# genome.  then write out all matches above pid threshold in that region
#
# good match
#   percent ID above pid_threshold
#   len above len_threshold
####################################################################################################
import sys


def geint(value, threshold):
    return int(value) >= threshold


def gefloat(value, threshold):
    return float(value) >= threshold


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    infile = sys.argv[1]
    sys.stderr.write('Blast search: {}\n'.format(infile))
    try:
        blast = open(infile, 'r')
    except OSError:
        sys.stderr.write('Cannot open blast file ({})'.format(infile))
        exit(1)

    # tests = {'length': {'column': 9, 'threshold': 20, 'fxn': geint},
    #          'pid': {'column': 10, 'threshold': 1.0, 'fxn': gefloat}}
    tests = {'length': {'column': 9, 'threshold': 100, 'fxn': geint},
             'pid': {'column': 10, 'threshold': 98.0, 'fxn': gefloat}}

    # subject(hit) begin and end colmns in tabular output
    query_id = 0
    subj_id = 4
    subj_begin = 6
    subj_end = 7
    al_len = 9

    # first pass, find the region for each transcript

    region = {}
    region_n = 0

    for line in blast:
        if line.startswith('#'):
            continue

        field = line.split("\t")

        status = 'success'
        for test in tests:
            thistest = tests[test]

            if thistest['fxn'](field[thistest['column']], thistest['threshold']):
                continue
            else:
                status = 'fail'
                break

        if status is 'success':

            qid = field[query_id]
            sid = field[subj_id]
            length = int(field[al_len])
            begin = int(field[subj_begin])
            end = int(field[subj_end])

            if begin > end:
                # make sure begin < end
                begin, end = end, begin

            if qid not in region:
                region[qid] = {}

            if sid in region[qid]:
                region[qid][sid]['begin'] = min(region[qid][sid]['begin'], begin)
                region[qid][sid]['end'] = max(region[qid][sid]['end'], end)
                region[qid][sid]['length'] = max(region[qid][sid]['length'], length)

            else:
                region[qid][sid] = {'begin': begin, 'end': end, 'length': length}
                region_n += 1

    # filter the regions and find those that have more than a minimum length
    sys.stdout.write('\nRegions ({}):\n'.format(region_n))
    selected = []
    selidx = {}
    minlen = 150
    segment = 100
    select_n = 0
    for q in region:
        # sys.stdout.write('{}\n'.format(q))

        for s in region[q]:
            begin = region[q][s]['begin']
            end = region[q][s]['end']
            regionlen = end - begin + 1
            sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format(q, s, begin, end, regionlen))

            if regionlen >= minlen and region[q][s]['length'] > segment:

                selected.append({'qid': q, 'sid': s, 'begin': begin, 'end': end})
                if q not in selidx:
                    selidx[q] = {}
                selidx[q][s] = selected[-1]
                select_n += 1

    sys.stdout.write('\nSelected Regions ({}):\n'.format(select_n))
    for sel in selected:
        sys.stdout.write('\t{}\t{}\t{}\t{}\t{}\n'.format(
            sel['qid'], sel['sid'], sel['begin'], sel['end'], sel['end'] - sel['begin'] + 1))

    sys.stdout.write('\nOverlapping Regions:\n')
    order = sorted(selected, key=lambda x: (x['sid'], x['begin']))
    pre = order[0]
    block = [pre]
    block_n = 0
    new = False

    for i in range(1, len(order)):
        cur = order[i]

        if cur['sid'] == pre['sid']:
            if cur['begin'] <= pre['end']:
                # sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                #     pre['qid'], pre['begin'], pre['end'],
                #     cur['qid'], cur['begin'], cur['end'], cur['sid']))
                block.append(cur)
            else:
                new = True
        else:
            new = True

        if new:
            # write out block
            sys.stdout.write('\nblock_{}\n'.format(block_n))
            for b in block:
                sys.stdout.write('\t{}\t{}\t{}\t{}\n'.format(
                    b['qid'], b['begin'], b['end'], b['sid']))
            block_n += 1
            block = [cur]
            new = False

        pre = cur

sys.stdout.write('\nblock_{}\n'.format(block_n))
for b in block:
    sys.stdout.write('\t{}\t{}\t{}\t{}\n'.format(
        b['qid'], b['begin'], b['end'], b['sid']))

# rewind file and process again pulling out regions of interest
sys.stdout.write('\nMapped transcripts:\n')
blast.seek(0, 0)
save = []
minlen = 100
buffer = 1000
for line in blast:
    if line.startswith('#'):
        continue

    field = line.split("\t")
    qid = field[query_id]
    sid = field[subj_id]
    if qid in selidx:
        if sid in selidx[qid]:
            region = selidx[qid][sid]
            begin = int(field[subj_begin])
            end = int(field[subj_end])
            if begin > end:
                begin, end = end, begin
                field[subj_begin] = str(begin)
                field[subj_end] = str(end)

            if (    (begin <= region['begin'] - buffer and end >= region['begin'] - buffer) or
                    (end >= region['end'] + buffer and begin <= region['end'] + buffer) or
                    (begin >= region['begin'] - buffer and end <= region['end'] + buffer)):
                save.append(field[:])

for hit in sorted(save, key=lambda x: (x[subj_id], int(x[subj_begin]))):
    # line = '\t'.join(hit)
    sys.stdout.write('{}'.format('\t'.join(hit)))

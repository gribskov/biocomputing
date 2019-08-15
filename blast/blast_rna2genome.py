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
    tests = {'length': {'column': 9, 'threshold': 30, 'fxn': geint},
             'pid': {'column': 10, 'threshold': 98.0, 'fxn': gefloat}}

    # subject(hit) begin and end colmns in tabular output
    query_id = 0
    subj_id = 4
    subj_begin = 6
    subj_end = 7

    # first pass, find the region for each transcipt

    region = {}

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

            # make sure begin < end
            begin = int(field[subj_begin])
            end = int(field[subj_end])
            if begin > end:
                begin, end = end, begin

            qid = field[query_id]
            sid = field[subj_id]

            if qid not in region:
                region[qid] = {}

            if sid in region[qid]:
                region[qid][sid]['begin'] = min(region[qid][sid]['begin'], begin)
                region[qid][sid]['end'] = max(region[qid][sid]['end'], end)

            else:
                region[qid][sid] = {'begin': begin, 'end': end, 'id': field[subj_id],
                                    'query': field[query_id]}

    for q in region:
        sys.stdout.write('{}\n'.format(q))

        for s in region[q]:
            begin = region[q][s]['begin']
            end = region[q][s]['end']

            sys.stdout.write('\t{}:{}\t{}\t{}\t{}\n'.
                             format(q, s, begin, end, end - begin + 1))

    # may be good to filter for regions over a minimum length

# rewind file and process again pulling out regions of interest
blast.seek(0, 0)
save = []
minlen = 100
for line in blast:
    if line.startswith('#'):
        continue

    field = line.split("\t")
    qid = field[query_id]
    sid = field[subj_id]
    if qid in region:
        if sid in region[qid]:
            if minlen > (region[qid][sid]['end'] - region[qid][sid]['begin']):
                print('skipping {}:{}\t{} - {} = {} < {}'.format(qid,sid,region[qid][sid]['end'],
                                                          region[qid][sid]['begin'],
                                                          region[qid][sid]['end'] -
                                                          region[qid][sid]['begin'] + 1, minlen))
                continue

            begin = int(field[subj_begin])
            end = int(field[subj_end])
            if begin > end:
                begin, end = end, begin
                field[subj_begin] = str(begin)
                field[subj_end] = str(end)

            if ((begin <= region[qid][sid]['begin'] and end >= region[qid][sid]['begin']) or
                (end >= region[qid][sid]['end'] and begin <= region[qid][sid]['end']) or
                (begin >= region[qid][sid]['begin'] and end <= region[qid][sid]['end'])):
                save.append(field[:])

for hit in sorted(save, key=lambda x: (x[subj_id], int(x[subj_begin]))):
    # line = '\t'.join(hit)
    sys.stdout.write('{}'.format('\t'.join(hit)))

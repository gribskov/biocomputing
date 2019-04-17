'''=================================================================================================
Cluster Trinity isoforms based on comparison to Uniref

Search is assumed to be diamond (blast tabular), tab separated
TRINITY_DN88428_c0_g1_i1   1146   1106   705   A0A059DJS1_EUCGR   290   1   135   135   65.2   4.2e-45   A0A059DJS1_EUCGR

qname qlen qbegin qend
sname slen sbegin send
alignlen score evalue stitle
================================================================================================='''
import sys
from blast import Blast

# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    infile = sys.argv[1]
    sys.stderr.write('Blast search: {}\n'.format(infile))
    blast = Blast(file=sys.argv[1])

    format = 'qname qlen qbegin qend sname slen sbegin send alignlen score evalue stitle'
    nfields = blast.setFormat(format)

    record = []
    record_n = 0
    sidx = {}
    qidx = {}
    threshold = 1e-5

    # read the search and store all matches over a threshold
    while blast.next():
        # print(blast.line)
        # print('query:{}\tsubject:{}'.format(blast.qname, blast.sname))

        if float(blast.evalue) <= threshold:
            d = blast.toDict()
            record.append(d)

            if blast.sname in sidx:
                sidx[blast.sname].append(d)
            else:
                sidx[blast.sname] = [d]

            if blast.qname in qidx:
                qidx[blast.qname].append(d)
            else:
                qidx[blast.qname] = [d]

            record_n += 1
            if record_n > 1000:
                break

    sys.stderr.write('{} records read from blast file\n'.format(record_n))

    # the subjects are the reference sequences, sort by number of matches
    cluster = {}
    for subj in sorted(sidx, key=lambda x: len(sidx[x]), reverse=True):
        print(subj)
        for hit in sidx[subj]:
            # print('\t{}'.format(hit))


#     for subj in sidx:
#         r = record[sidx[subj][0]]['stitle']
#         l = record[sidx[subj][0]]['slen']
#         print('{}\tlen={}\t{}'.format(subj, l, r))
#         for i in sorted(sidx[subj], key=lambda x: record[x]['qname']):
#             q = record[i]
#             qcov = (q['qend'] - q['qbegin'] + 1) / q['qlen']
#             scov = (q['send'] - q['sbegin'] + 1) / q['slen']
#             print('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
#                 scov, q['sbegin'], q['send'], q['slen'],
#                 q['qname'], qcov, q['qbegin'], q['qend'], q['qlen'], q['evalue']))

exit(0)

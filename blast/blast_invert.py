'''=================================================================================================
Invert a blast search - reorganize by subject instead of query.  This is usefule for transcript
sequences to see which predicted transcripts correspond to the same gene

Search is assumed to be diamond (blast tabular)
TRINITY_DN88428_c0_g1_i1        1146    1106    705
A0A059DJS1_EUCGR        290     1       135
135     65.2    4.2e-45 A0A059DJS1_EUCGR

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
    sys.stderr.write('Blast serach: {}\n'.format(infile))
    blast = Blast(file=sys.argv[1])

    format = 'qname qlen qbegin qend sname slen sbegin send alignlen score evalue stitle'
    nfields = blast.setFormat(format)

    record = []
    record_n = 0
    sidx = {}
    while blast.next():
        # print(blast.line)
        # print('query:{}\tsubject:{}'.format(blast.qname, blast.sname))
        d = blast.toDict()

        if blast.sname in sidx:
            sidx[blast.sname].append(record_n)
        else:
            sidx[blast.sname] = [record_n]

        record.append(d)
        record_n += 1

    sys.stderr.write('{} records read from blast file\n'.format(record_n))

    for subj in sidx:
        print('{}'.format(subj))
        for i in sidx[subj]:
            q = record[i]
            print(
                '\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(q['qname'], q['sbegin'], q['send'], q['slen'],
                                                      q['qbegin'], q['qend'], q['qlen']))

exit(0)

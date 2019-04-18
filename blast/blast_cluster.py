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


class Cluster:
    """---------------------------------------------------------------------------------------------
    groups of pairs of sequences with blast edges
    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.seqlist = {}
        self.group = []

    def group_new(self, record):
        """-----------------------------------------------------------------------------------------
        Create a new group

        :param record: dict, blast record for a subj-query edge
        :return: int, size of group list
        -----------------------------------------------------------------------------------------"""
        group = len(self.group)
        self.group.append([record])
        self.seqlist[record['sname']] = group
        self.seqlist[record['qname']] = group

        return group

    def group_num(self, record):
        """-----------------------------------------------------------------------------------------
        check the subject and query and see if they are in a group

        :param record: dict, blast record for a subj-query edge
        :return: tuple of int, group indices
        -----------------------------------------------------------------------------------------"""
        sgroup = None
        if record['sname'] in self.seqlist:
            sgroup = self.seqlist[record['sname']]

        qgroup = None
        if record['qname'] in self.seqlist:
            qgroup = self.seqlist[record['qname']]

        return (sgroup, qgroup)

    def add(self, group, record):
        """-----------------------------------------------------------------------------------------
        Add a new query to an existing group

        :param sgroup: number of subject group
        :param record: dict, blast record for a subj-query edge
        :return: int, edges in group
        -----------------------------------------------------------------------------------------"""
        self.group[group].append(record)
        self.seqlist[record['qname']] = group
        self.seqlist[record['sname']] = group

        return len(self.group[group])

    def group_merge(self, keep, discard):
        """-----------------------------------------------------------------------------------------
        merge 2 existing groups.  The discard group is merged into the keep group

        :param keep:
        :param discard:
        :return: int, edges in merged group
        -----------------------------------------------------------------------------------------"""
        for edge in self.group[discard]:
            self.group[keep].append(edge)
            self.seqlist[edge['sname']] = keep
            self.seqlist[edge['qname']] = keep

        self.group[discard] = None
        return len(self.group[sgroup])


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
    cluster = Cluster()

    for subj in sorted(sidx, key=lambda x: len(sidx[x]), reverse=True):
        print(subj)
        for edge in sidx[subj]:
            (sgroup, qgroup) = cluster.group_num(edge)

            if sgroup is None:
                if qgroup is None:
                    # both unknown, start new group
                    sgroup = cluster.group_new(edge)
                    qgroup = sgroup
                    continue

                else:
                    # sgroup unknown, qgroup known, add s to q group
                    cluster.add(qgroup,edge)
                    continue

            if qgroup is None:
                # sgroup is known, qgroup is unknown
                cluster.add(sgroup, edge)
                qgroup = sgroup
                continue

            if sgroup != qgroup:
                # merge groups
                cluster.group_merge(sgroup, qgroup)

    for group in cluster.group:
        if group is None:
            continue

        print(group[0]['qname'])
        for edge in group:
            print("\t", edge['stitle'])
exit(0)

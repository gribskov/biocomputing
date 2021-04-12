"""=================================================================================================


Michael Gribskov     12 April 2021
================================================================================================="""
import sys


class SingleLinkage:
    """=============================================================================================

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Data should be entered as a dictionary one item at a time. Each item of distance data is
        stored as a list with the fields determined by the label and keys.
            labelcol - list of keys in input data to use as labels (typically two fields)
            keys - list of keys in input data to include in data list
            kidx - conversion from key bact to original column name
            data = [label1, label2, key1, key2, key3, ...]

        Data is stored in the barrio dictionary as
            barrio[label1] = [[data1], [data2] ...]

        -----------------------------------------------------------------------------------------"""
        self.barrio = {}
        self.group = []
        self.groupidx = {}
        self.names = {}
        self.labels = {}
        self.keys = {}
        self.kidx = []

    def append(self, datadict):
        """-----------------------------------------------------------------------------------------
        create the data entry list and store in barrio under labels in labelcol

        :param datadict: dict with distance data
        :return:
        -----------------------------------------------------------------------------------------"""
        data = []
        for k in self.keys:
            if k in self.labels:
                name = datadict[k]
                if name not in self.names:
                    self.names[name] = len(self.names)
                    self.barrio[self.names[name]] = []
                data.append(self.names[name])
            else:
                data.append(datadict[k])

        for l in self.labels:
            name = self.names[datadict[l]]
            self.barrio[name].append(data)

        return None

    @staticmethod
    def invert_names(names):
        """-----------------------------------------------------------------------------------------
        return a dictionary with the keys and values reverse so that numerical IDs can be
        converted back to strings

        :param names: dict, usually self.names
        :return: dict
        -----------------------------------------------------------------------------------------"""
        reverse = {}
        for id in names:
            reverse[names[id]] = id

        return reverse

    def set_keys(self, keys):
        """-----------------------------------------------------------------------------------------
        Set up the coversion between data column names and indices.  Include labels in the keys

        :param keys: list, column names in original data
        :return: int, number keys in kidx
        -----------------------------------------------------------------------------------------"""
        for k in keys:
            if k in self.keys:
                # known key, skip
                continue

            self.keys[k] = len(self.kidx)
            self.kidx.append(k)

        return len(self.kidx)

    def single(self):
        """-----------------------------------------------------------------------------------------
        second try using the structure of the data
        :return:
        -----------------------------------------------------------------------------------------"""
        gidx = self.groupidx
        for l1 in self.barrio:
            known = set()
            unknown = set()
            if l1 in gidx:
                known.add(l1)
            else:
                unknown.add(l1)

            for b in self.barrio[l1]:
                l2 = b[1]
                if l2 in gidx:
                    known.add(l2)
                else:
                    unknown.add(l2)

            # print('{}:{}:{}'.format(l1, known, unknown))
            if not known:
                # all are unknown, make a new group and assign all to it
                g = len(self.group)
                self.group.append([])
            else:
                g = self.groupidx[known.pop()]

            for member in unknown:
                self.group[g].append(member)
                gidx[member] = g
            for member in known:
                if self.groupidx[member] != g:
                    self.group.g += self.group[groupidx[member]]
                    gidx.member = g

            # print('{}:{}'.format(g, self.group[g]))

        return

    def write_groups(self, fh, names=True):
        """-----------------------------------------------------------------------------------------
        write all groups to filehandle

        :return:
        -----------------------------------------------------------------------------------------"""
        if names:
            nidx = SingleLinkage.invert_names(self.names)

            for g in range(len(self.group)):
                fh.write('group {}: '.format(g))
                members = []
                for each in self.group[g]:
                    members += nidx[each]
                fh.write('{}\n'.format(', '.join(members)))

        else:
            for g in range(len(self.group)):
                fh.write('{}: {}\n'.format(g, self.group[g]))

        return


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    cluster = SingleLinkage()
    cluster.set_keys(['sname', 'qname', 'evalue', 'id'])
    cluster.labels = ['sname', 'qname']

    d0 = [{'sname':'a', 'qname':'b', 'evalue':1e-100, 'id':80},
          {'sname':'a', 'qname':'c', 'evalue':1e-100, 'id':80},
          {'sname':'b', 'qname':'c', 'evalue':1e-100, 'id':80},
          {'sname':'d', 'qname':'e', 'evalue':1e-100, 'id':80},
          ]

    for d in d0:
        cluster.append(d)

    cluster.single()

    d1 = [{'sname':'a', 'qname':'b', 'evalue':1e-100, 'id':80},
          {'sname':'a', 'qname':'c', 'evalue':1e-100, 'id':80},
          {'sname':'b', 'qname':'c', 'evalue':1e-100, 'id':80},
          {'sname':'d', 'qname':'e', 'evalue':1e-100, 'id':80},
          {'sname':'f', 'qname':'g', 'evalue':1e-100, 'id':80},
          {'sname':'f', 'qname':'b', 'evalue':1e-100, 'id':80},
          {'sname':'d', 'qname':'a', 'evalue':1e-100, 'id':80},
          ]

    cluster = SingleLinkage()
    cluster.set_keys(['sname', 'qname', 'evalue', 'id'])
    cluster.labels = ['sname', 'qname']
    for d in d1:
        cluster.append(d)

    cluster.single()
    cluster.write_groups(sys.stdout)

    d1 = [{'sname':'a', 'qname':'b', 'evalue':1e-50, 'id':80},
          {'sname':'a', 'qname':'c', 'evalue':1e-100, 'id':80},
          {'sname':'b', 'qname':'c', 'evalue':1e-100, 'id':80},
          # {'sname':'b', 'qname':'h', 'evalue':1e-100, 'id':80},
          {'sname':'c', 'qname':'d', 'evalue':1e-100, 'id':80},
          {'sname':'c', 'qname':'e', 'evalue':1e-100, 'id':80},
          {'sname':'d', 'qname':'e', 'evalue':1e-100, 'id':80},
          {'sname':'e', 'qname':'f', 'evalue':1e-50, 'id':80},
          # {'sname':'f', 'qname':'g', 'evalue':1e-100, 'id':80},
          {'sname':'g', 'qname':'h', 'evalue':1e-100, 'id':80},
          {'sname':'g', 'qname':'j', 'evalue':1e-100, 'id':80},
          {'sname':'g', 'qname':'l', 'evalue':1e-100, 'id':80},
          {'sname':'h', 'qname':'i', 'evalue':1e-100, 'id':80},
          {'sname':'i', 'qname':'j', 'evalue':1e-100, 'id':80},
          {'sname':'j', 'qname':'k', 'evalue':1e-100, 'id':80},
          {'sname':'k', 'qname':'l', 'evalue':1e-100, 'id':80},
          ]

    cluster = SingleLinkage()
    cluster.set_keys(['sname', 'qname', 'evalue', 'id'])
    cluster.labels = ['sname', 'qname']
    for d in d1:
        cluster.append(d)

    cluster.single()
    cluster.write_groups(sys.stdout)

    exit(0)

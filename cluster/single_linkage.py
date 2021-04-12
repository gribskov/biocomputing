"""=================================================================================================


Michael Gribskov     12 April 2021
================================================================================================="""


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
        self.groups = []
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

    exit(0)

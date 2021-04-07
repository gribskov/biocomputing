class box(object):

    def __init__(self):
        self.content = None
        self.i = None

    def __iter__(self):
        for item in self.content:
            yield item

    def __len__(self):
        return len(self.content)

    def __repr__(self):
        return super().__repr__()

    def __str__(self):
        return ','.join(self.content)

if __name__ == '__main__':

    caja = box()
    caja.content = ['doll', 'truck']
    caja.content.append('lightsaber')

    print(caja)
    print(caja.__repr__())
    print(caja.__str__())

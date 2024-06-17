"""=================================================================================================


Michael Gribskov     12 February 2024
================================================================================================="""
from sys import getsizeof

class SequenceIterator:
    def __init__(self, sequence):
        self._sequence = sequence
        self._index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._index < len(self._sequence):
            item = self._sequence[self._index]
            self._index += 1
            return item
        else:
            raise StopIteration

# for i in range(99, 102):
#     if i >= 100:
#         raise Exception('TooBigError')



a = [1,3,5]
iter_a = SequenceIterator(a)

for i in (0, 1):
    print(next(iter_a))

print(next(iter_a))
# print(next(iter_a))

million = [i for i in range(1000000)]
for x in million:
    pass

print(f'million is {len(million)} long')

class SequenceIterator:
    def __init__(self, limit):
        self._limit = limit
        self._index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._index < self._limit:
            item = self._index
            self._index += 1
            return item
        else:
            raise StopIteration

million = SequenceIterator(1000000)
b = iter(million)
for x in b:
    i = i + 1

print(f'million is {getsizeof(b)} long')

exit(0)


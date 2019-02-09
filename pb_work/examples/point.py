"""=================================================================================================
Point class example

Michael Gribskov     09 February 2019
================================================================================================="""
import math


class Point:

    def __init__(self, x, y):
        self._pos = [x, y]

    def magnitude(self):
        return math.sqrt(self._pos[0] ** 2 + self._pos[1] ** 2)

    def add(self, other):
        self._pos[0] += other._pos[0]
        self._pos[1] += other._pos[1]

        return self.magnitude()

    def values(self):
        return self._pos


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    pt1 = Point(1, 2)
    pt2 = Point(4, 4)

    print('pt1 length = {}'.format(pt1.magnitude()))
    print('pt2 length = {}'.format(pt2.magnitude()))

    length = pt1.add(pt2)
    print('pt1 + pt2 = {}'.format(length))

    x, y = pt1.values()
    print('pt1 length = {}\tx = {}\ty = {}'.format(pt1.magnitude(), x, y))

    xy = []
    xy = pt1.values()
    print('pt1 length = {}\tx = {}\ty = {}'.format(pt1.magnitude(), xy[0], xy[1]))
    xy[0] = 10
    xy = pt1.values()
    print('pt1 length = {}\tx = {}\ty = {}'.format(pt1.magnitude(), xy[0], xy[1]))

    exit(0)

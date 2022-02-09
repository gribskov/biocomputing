"""=================================================================================================
Point class example

Michael Gribskov     09 February 2019
================================================================================================="""
import math


class Point:

    def __init__(self, x=0, y=0):
        self._pos = []
        self.setvalues(x,y)

    def magnitude(self):
        return math.sqrt(self._pos[0] ** 2 + self._pos[1] ** 2)

        self.setvalues(x, y)

    def setvalues(self, x, y):
        self._pos = [x, y]

        return True

    def getvalues(self):
        return self._pos[:]

    def add(self, other):
        self._pos[0] += other._pos[0]
        self._pos[1] += other._pos[1]

        return self.magnitude()

    def add_better(self, other):
        x0, y0 = self.getvalues()
        x1, y1 = other.getvalues()
        self._pos = [x0 + x1, y0 + y1]

        return self.magnitude()

class Shape(Point):

    def __init__(self, x=0, y=0, shape='circle', size=1):
        """
        Shape subclass adds shape and size to point
        for now shapes are limited to circle and square
        size can be in any arbitrary plotting units compatible with the plot
        """
        self.shape = shape
        self.size = size
        self = super().__init__(x,y)



class Polar(Point):
    """ _pos = [angle, radius]"""

    # Note: __init__ inherited from Point

    def magnitude(self):
        return self._pos[1]

    def to_cartesian(self):
        return self._pos[1] * math.cos(math.radians(self._pos[0])), \
                self._pos[1] * math.sin(math.radians(self._pos[0]))


    def add(self, other):
        x1, y1 = self.to_cartesian()
        x2, y2 = other.to_cartesian()

        self._pos[0] = math.degrees(math.atan2(y1 + y2, x1 + x2))
        self._pos[1] = math.sqrt((x1 + x2) ** 2 + (y1 + y2) ** 2)

        return self._pos[1]


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # po1 = Polar(0.0, 1.0)
    # print('po1 values = {}'.format(po1.values()))
    # po2 = Polar(90.0, 1.0)
    # print('po2 values = {}'.format(po2.values()))
    # po1.add(po2)
    # print('po1 + po2 = {}'.format(po1.values()))

    pt1 = Point(1, 2)
    pt2 = Point(4, 4)

    print('pt1 length = {}'.format(pt1.magnitude()))
    print('pt2 length = {}'.format(pt2.magnitude()))

    length = pt1.add(pt2)
    print('pt1 + pt2 = {}'.format(length))

    x, y = pt1.getvalues()
    print('pt1 length = {}\tx = {}\ty = {}'.format(pt1.magnitude(), x, y))

    xy = []
    xy = pt1.getvalues()
    print('pt1 length = {}\tx = {}\ty = {}'.format(pt1.magnitude(), xy[0], xy[1]))


    print('\nchanging xy, pt1 should not change')
    xy[0] = 10
    xy = pt1.getvalues()
    print('pt1 length = {}\tx = {}\ty = {}'.format(pt1.magnitude(), xy[0], xy[1]))

    s = Shape()
    s = Shape(1,2, 'square')
    print(f'object:{s}   shape={s.shape}   size={s.size}   pos={s.getvalues()}   length={s.magnitude():.3f}')

    exit(0)

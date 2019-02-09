"""=================================================================================================
Point class example

Michael Gribskov     09 February 2019
================================================================================================="""
import math

class Point:

    def __init__(self):
        self.pos = []

    def magnitude(self):
        return math.sqrt(self.pos[0] ** 2 + self.pos[1] ** 2)


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    pt1 = Point()
    pt1.pos = [1, 2]
    pt2 = Point()
    pt2.pos = [4,4]

    print('pt1 length = {}'.format(pt1.magnitude()))
    print('pt2 length = {}'.format(pt2.magnitude()))


    exit(0)

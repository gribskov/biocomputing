"""=================================================================================================
from https://realpython.com/primer-on-python-decorators/

Michael Gribskov     30 March 2021
================================================================================================="""


class Circle:

    def __init__(self, radius):
        self._radius = radius

    @property
    def radius(self):
        """Get value of radius"""
        r = '{:.3f}'.format(self._radius)
        return float(r)

    @radius.setter
    def radius(self, value):
        """Set radius, raise error if negative"""
        if value >= 0:
            self._radius = value
        else:
            raise ValueError("Radius must be positive")

    @property
    def area(self):
        """Calculate area inside circle"""
        return self.pi() * self.radius ** 2

    def cylinder_volume(self, height):
        """Calculate volume of cylinder with circle as base"""
        return self.area * height

    @classmethod
    def unit_circle(cls):
        """Factory method creating a circle with radius 1"""
        return cls(1)

    @staticmethod
    def pi():
        """Value of π, could use math.pi instead"""
        return 3.1416


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    c = Circle(5.1234567)

    print('radius', c.radius )
    print('_radius', c._radius)
    try:
        c.radius = -1
    except ValueError:
        print('negative radius')
        pass
    print('-radius', c.radius)  # negative radius, should leave radius unchanged
    print( 'area', c.area)

    c.radius = 2
    print( 'area', c.area)

    print('cylinder', c.cylinder_volume(height=4))

    print(c)
    c = Circle.unit_circle()
    print(c)
    print('radius', c.radius)

    print('pi', Circle.pi(), c.pi())

    # try to set area, should fail
    c.area = 100


    exit(0)

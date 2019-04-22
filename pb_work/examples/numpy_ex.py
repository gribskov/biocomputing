import numpy
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    dist = [[10, 10, 2],
            [5, 5, 1],
            [10, 5, 2]]
    n = 1000

    data = []
    for center in dist:
        data += numpy.random.normal(loc=[center[0], center[1]], scale=center[2],
                                    size=[n, 2]).tolist()

    print('len', len(data))
    x = []
    y = []
    for point in data:
        x.append(point[0])
        y.append(point[1])

    plt.scatter(x, y, marker='.', linewidth=1.0)
    plt.show()

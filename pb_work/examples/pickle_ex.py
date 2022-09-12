"""=================================================================================================
Pickle examples

Michael Gribskov     19 April 2021
================================================================================================="""
import pickle

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    dummy = {1:'A', 2:'B', 3:'C', 4:'D', 5:'E'}
    topickle = open('test.pickle', 'wb')
    pickle.dump(dummy, topickle)
    topickle.close()

    frompickle = open('test.pickle', 'rb')
    dummy = pickle.load(frompickle)
    print(dummy)

    dummy = {1:'A', 2:'B', 3:'C', 4:'D', 5:'E'}
    dummy[6] = [{'name':'level1', 'data':[5, 4, 3, 2, 1]},
                {'name':'level1', 'data':[5, 4, 3, 2, 1]}]
    topickle = open('test.pickle', 'wb')
    pickle.dump(dummy, topickle)
    topickle.close()

    frompickle = open('test.pickle', 'rb')
    dummy = pickle.load(frompickle)
    print(dummy)

    # circular
    circ = {1:'A', 2:'B', 3:'C', 4:'D'}
    circ[5] = circ
    print(f'\ncirc:{circ}\n nested circ:{circ[5]}')
    topickle = open('test.pickle', 'wb')
    pickle.dump(circ, topickle)
    topickle.close()

    frompickle = open('test.pickle', 'rb')
    circ = pickle.load(frompickle)
    print(f'circ:{circ}\n nested circ:{circ[5]}')

    # object
    class data:
        def __init__(self):
            self.a = None
            self.b = None
            self.c = None

    test = data()
    test.a = [1,2,3,4,5]
    test.b = {'a':100, 'b':1000}
    test.c = data()
    print(f'test:{test}\ntest.a:{test.a}  test.b:{test.b}  test.c:{test.c}')
    topickle = open('test.pickle', 'wb')
    pickle.dump(test, topickle)
    topickle.close()

    frompickle = open('test.pickle', 'rb')
    test = pickle.load(frompickle)
    print(f'test:{test}\ntest.a:{test.a}  test.b:{test.b}  test.c:{test.c}')
    frompickle.close()

    exit(0)

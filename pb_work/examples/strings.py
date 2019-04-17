"""=================================================================================================
String examples

Michael Gribskov     14 April 2019
================================================================================================="""

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    name = 'Bob'
    age = 23

    print('name = {}\tage = {}'.format(name, age))
    print('name = {0}\tage = {1}'.format(name, age))
    print('age = {1}\tname = {0}'.format(name, age))
    print('name = {name}\tage = {age}'.format(name=name, age=age))

    bob = ['Bob', '23']
    ted = ['Ted', '21']
    carol = ['Carol', '24']

    print('\nname = {0}\tage = {1}'.format(ted[0], ted[1]))
    print('name = {0[0]}\tage = {0[1]}'.format(bob))
    print('name = {0[0]}\tage = {1[1]} ({1[0]})'.format(bob, carol))

    person = {'name': 'Bob', 'age': 23}
    print('\nname = {name}\tage = {age}'.format(name=person['name'], age=person['age']))
    print('name = {name}\tage = {age}'.format(**person))

    print('\nf-strings')
    name = 'Bob'
    age = 23
    print(f'name = {name}\tage = {age}')
    print(f'name = {str.lower(name)}\tage = {age}')

    exit(0)

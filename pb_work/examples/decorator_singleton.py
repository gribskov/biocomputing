"""=================================================================================================


Michael Gribskov     02 June 2021
================================================================================================="""
import functools

def singleton(cls):
    """Make a class a Singleton class (only one instance)"""
    @functools.wraps(cls)
    def wrapper_singleton(*args, **kwargs):
        if not wrapper_singleton.instance:
            wrapper_singleton.instance = cls(*args, **kwargs)
        return wrapper_singleton.instance
    wrapper_singleton.instance = None
    return wrapper_singleton

@singleton
class TheOne:
    pass

@singleton
class TheTwo:
    pass

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    first_one = TheOne()
    another_one = TheOne()
    the_two = TheTwo()
    the_tutu = TheTwo()
    third_one = TheOne()

    print(id(first_one))
    print(id(another_one))
    print(id(third_one))
    print(id(the_two))
    print(id(the_tutu))

    exit(0)

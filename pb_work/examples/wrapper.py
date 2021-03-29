"""=================================================================================================


Michael Gribskov     29 March 2021
================================================================================================="""


def wrapper(func):
    def inner():
        print("Something before the function is called.")
        func()
        print("Something after the function is called.\n")

    return inner


def say_whee():
    print("Whee!")


def my_decorator(func):
    def inner():
        print("Something before the function is called.")
        func()
        print("Something after the function is called.\n")

    return inner


@my_decorator
def deco_whee():
    print('deco whee')


def do_twice(func):
    def wrapper_do_twice():
        func()
        func()

    return wrapper_do_twice


@do_twice
def greet(name):
    print(f'Hello {name}')

def do_twice_var(func):
    def wrapper_do_twice(*args, **kwargs):
        func(*args, **kwargs)
        func(*args, **kwargs)
    return wrapper_do_twice

@do_twice_var
def return_greeting(name):
    print('Creating greeting')
    return f'Hi {name}'

def do_twice_ret(func):
    def wrapper_do_twice(*args, **kwargs):
        func(*args, **kwargs)
        val = func(*args, **kwargs)
        return val
    return wrapper_do_twice

@do_twice_ret
def return_greeting(name):
    print('Creating greeting')
    return f'Hi {name}'
#--------------------------------------------------
# main
#--------------------------------------------------
if __name__ == '__main__':
    print('\nwrapper')
    whee = wrapper(say_whee)
    print(whee)
    whee()

    print('\ndecorator')
    print(deco_whee)
    deco_whee()

    # greet('World')

    hi_adam = return_greeting('Adam')
    print(hi_adam)



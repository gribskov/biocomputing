"""=================================================================================================
Some examples of decorators

Michael Gribskov     14 April 2019
================================================================================================="""


def hello(name):
    return f'Hello, {name}.'


def hola(name):
    return f'Hola, {name}.'


def greet(greeting_function, name):
    return greeting_function(name)


def language(idioma):
    if idioma == 'english':
        return hello
    elif idioma == 'spanish':
        return hola
    else:
        print(f'unknown language ({idioma})')
        exit(1)


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


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    print(hello('Bob'))

    name = 'Eric'
    print(hola('Eric'))

    print('\ngreeting function')
    for name in ('Bob', 'Eric'):
        print(greet(hello, name))
        print(greet(hola, name))

    print('\nusing returned function')
    for tongue in ('english', 'spanish'):
        greeting = language(tongue)
        print(greet(greeting, 'Bob'))

    print('\nwrapper')
    whee = wrapper(say_whee)
    print(whee)
    whee()

    print('\ndecorator')
    print(deco_whee)
    deco_whee()


    exit(0)

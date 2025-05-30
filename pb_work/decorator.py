"""=================================================================================================


Michael Gribskov     10 April 2025
================================================================================================="""

import time

def do_twice(func):
    def wrapper_do_twice(*args, **kwargs):
        func(*args, **kwargs)
        func(*args, **kwargs)

    return wrapper_do_twice


@do_twice
def return_greeting(name):
    print('Creating greeting')
    return f'Hi {name}'



def timer(func):
    """Print the runtime of the decorated function"""

    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()  # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()    # 2
        run_time = end_time - start_time
        print(f'Finished {func.__name__!r} in {run_time:.4f} secs')
        return value

    return wrapper_timer

@timer
def square_for(num_times):
    """simple sum of squares numtimes times with for loop"""
    s = 0
    for i in range(num_times):
        s += i ** 2
    return s

@timer
def square_while(num_times):
    """simple sum of squares numtimes times with while loop"""
    s = 0
    i = 0
    while i < num_times:
        s += i ** 2
        i += 1
    return s

@timer
def square_comp(num_times):
    """sum of squares with comprehension"""
    return sum([i ** 2 for i in range(num_times)])


def sqsub(i):
    """square function"""
    return i ** 2

def sqsubtimes(i):
    """square by multiplyng by self"""
    return i * i

@timer
def square_sqsub(num_times):
    """just like square_for with subroutine"""
    s = 0
    for i in range(num_times):
        s += sqsub(i)
    return s

@timer
def square_sqsubtimes(num_times):
    s = 0
    for i in range(num_times):
        s += sqsubtimes(i)
    return s

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # hi_adam = return_greeting('Adam')
    # print(hi_adam)


    i = 10000000
    square_for(i)
    square_while(i)
    square_comp(i)

    # sqsub(i)
    # sqsubtimes(i)
    square_sqsub(i)
    square_sqsubtimes(i)

    exit(0)

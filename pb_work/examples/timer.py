"""=================================================================================================
Use a decorator to time functions

Michael Gribskov     29 March 2021
================================================================================================="""
import time


def timer(func):
    """Print the runtime of the decorated function"""

    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()  # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_timer


@timer
def square_for(numtimes):
    """simple sum of squares numtimes times with for loop"""
    s = 0
    for i in range(numtimes):
        s += i ** 2
    return s


@timer
def square_while(num_times):
    """simple sum of squares numtimes times with while loop"""
    s = 0
    i = 0
    while i < numtimes:
        s += i ** 2
        i += 1
    return s


@timer
def square_comp(num_times):
    """sum of squares with comprehension"""
    return sum([i ** 2 for i in range(numtimes)])


def sqsub(i):
    """square function"""
    return i ** 2


def sqsubtimes(i):
    """square by multiplyng by self"""
    return i * i


@timer
def square_sqsub(num_times):
    s = 0
    for i in range(numtimes):
        s += sqsub(i)
    return s


@timer
def square_sqsubtimes(num_times):
    s = 0
    for i in range(numtimes):
        s += sqsubtimes(i)
    return s


if __name__ == '__main__':
    numtimes = 10000000
    square_for(numtimes)
    square_while(numtimes)
    square_comp(numtimes)
    square_sqsub(numtimes)
    square_sqsubtimes(numtimes)

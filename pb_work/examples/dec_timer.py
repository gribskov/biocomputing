"""=================================================================================================
A decorator that implements a timer

Michael Gribskov     14 April 2019
================================================================================================="""
import time


def timer(func):
    """Print the runtime of the decorated function"""
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value
    return wrapper_timer

@timer
def square_comp(num_times):
    s = 0
    for _ in range(num_times):
        s = sum([i**2 for i in range(10000)])

    return sum

@timer
def square_simple(num_times):
    s = 0
    for _ in range(num_times):
        for i in range(10000):
            s += i **2

    return s

@timer
def square_simple2(num_times):
    s = 0
    for i in range(10000):
        for _ in range(num_times):
            s += i ** 2

    return s

def sqsub(i,num_times):
    s = 0
    for _ in range(num_times):
        s +=  i **2

    return s

@timer
def square_sub(num_times):
    s = 0
    for i in range(10000):
        s += sqsub(i,num_times)

    return s

def sqsub2(i,num_times):
       return num_times * i ** 2

@timer
def square_sub2(num_times):
    s = 0
    for i in range(10000):
        s += sqsub2(i,num_times)

    return s
# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    square_simple(10000)
    square_comp(10000)
    square_simple2(10000)
    square_sub(10000)
    square_sub2(10000)

    exit(0)
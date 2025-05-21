"""=================================================================================================


Michael Gribskov     21 May 2025
================================================================================================="""
import sys
from time import perf_counter as time
from collections import defaultdict

class Primes:

    @staticmethod
    def stream():
        stack = defaultdict(lambda: [])
        yield 2
        yield 3
        stack[6] = [3]
        i = 5
        while True:
            next_composite = min(stack)
            while i < next_composite:
                # add primes up to the next composite
                yield i
                # stack[i * i].append(i * 2)
                stack[i * 2].append(i)
                i += 2


            # update stack with composite
            for val in stack[next_composite]:
                stack[val + next_composite].append(val)
            del stack[next_composite]
            if i == next_composite:
                i += 2

            # sys.stderr.write(f'{i}  min={min(stack)}\n')

        return


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
def verify(from_n, *vals):
    stream = Primes.stream()
    for _ in range(from_n): next(stream)
    ok = True
    for v in vals:
        ns = next(stream)
        if v != ns:
            ok = False
            break

    # print(f'test:{from_n}:{','.join(str(v) for v in vals)}\t{ok}')

elapse = 0.0
elapse_min = 10000
elapse_max = 0
repeat = 20
for i in range(repeat):
    start = time()
    verify(0, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
    verify(10, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71)
    verify(100, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601)
    verify(1000, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017)
    end = time()
    delta = end - start
    elapse += end - start
    elapse_min = min(elapse_min, delta)
    elapse_max = max(elapse_max, delta)

print(f'min:{1000*elapse_min:.3f}  ave:{1000*(elapse/repeat):.3f}  max:{1000*elapse_max:.3f}')

exit(0)

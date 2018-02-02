# from the sololearn tutorial
def add_five(x):
    return x + 5


nums = [11, 22, 33, 44, 55]
result = list(map(add_five, nums))
print(result)

# We could have achieved the same result more easily by using lambda syntax.
# they say

result = list(map(lambda x: x + 5, nums))
print(result)


# i say, we could do it better with a properly designed function

def add_five_list(x):
    for i in range(len(x)):
        x[i] += 5
    return x

result = add_five_list(nums)
print(result)


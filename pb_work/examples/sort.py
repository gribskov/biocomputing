words = "This is a test string from Andrew".split()
print(words)
print(sorted(words, key=str.lower))  # case insensitive
print(sorted(words, key=lambda word: word[-1]))  # sort by last letter
print(sorted(words, key=lambda word: word[::-1]))  # sort in reverse order

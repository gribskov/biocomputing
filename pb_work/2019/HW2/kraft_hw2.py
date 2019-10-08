# loads fastq file
try:
    infile = open('../HW1/8044.5k.fastq', 'r')
except:
    print('error loading file')
# makes fastq file into a list, creates empty lists base_string and q_string where base call and Q score will be
# stored after stripping '\n' from list
seq = list(infile)
base_list = []
q_list = []

# stores Q score in q_list after stripping white space
for i in seq[3::4]:
    q_list.append(i.strip())

# stores base call in base_list after stripping white space
for c in seq[1::4]:
    base_list.append(c.strip())

# turns base_list and q_list from a list to string
base_string = ''.join(base_list)
q_string = ''.join(q_list)

# counts number of sequences and number of bases
n_line = 0
length = 0
num_base = 0
for line in seq:
    n_line += 1
    if n_line % 4 == 2:
        length += 1

print('Number of sequences =', length, '\n', 'Number of bases =', len(base_string))
# function to convert the quality character string to a list of quality integer values
def quality(quality_string):
    for i in quality_string:
        return ord(i) - 33

"""takes ASCII character value and uses this to find the integer value. Adding 33 to integer value gives the
Q value. Example: print (quality('A') prints out 33 """
# Creates an empty list to store quality score output
q_output = []

# Iterate through quality score string to find corresponding Q value, stores values in the empty list q_ouput
for i in q_string:
    q_output.append([quality([i])])
base = []
for i in base_string:
        base.extend([i])
        # function to calculate A,C,G,T content of each sequence

def composition(sequence_string):
    bases = dict()
    base_total = 0
    for i in sequence_string:
        if i in bases:
            bases[i] += 1
        else:
            bases[i] = 1
    for b in sequence_string:
        base_total += 1
    print('The number of bases are:', (bases), 'The total number of bases are:', (base_total))
"""function used to calculate the A,C,G,T content for the sequences"""
composition(base_string)  # runs composition
# Counts number of each base type in sequence, used to check output of composition function
# A = 0
# C = 0
# G = 0
# T = 0
# N = 0
#
# for i in base_string:
#     if i == 'A':
#         A += 1
#     if i == 'C':
#         C += 1
#     if i == 'G':
#         G += 1
#     if i == 'T':
#         T += 1
#     if i == 'N':
#         N += 1
#
# print('The number of different bases are as follows:', '\n''A =', A, 'C =', C, 'G =', G, 'T =', T, 'N =', N)
# function to truncate sequence at first base with a quality score lower than a threshold, our threshold is 20
trimmed = []
for i in q_output:
    if [i] in [q_output] >= 20:
        trimmed.extend([base_string[i]])
def fastqNext(fh):
    n = 0
    seq = ''

    for line in fh:
        n += 1
        line = line.rstrip()
        if n % 4 == 2:
            seq = line
        if n == 4:
            break
    return seq


fh = open('small.fq', 'r')

base = ['A', 'C', 'G', 'T']

# n = 0
# nbase = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
# for line in fh:
#     n += 1
#     line = line.rstrip()
#     if n % 4 == 2:
#         print(line)
#         for b in line:
#             nbase[b] += 1
#
#         for b in nbase.keys():
#             print(b, nbase[b])

while True:
    seq = fastqNext(fh)
    print(seq)

    if not seq:
        break

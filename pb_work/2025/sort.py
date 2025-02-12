"""=================================================================================================


Michael Gribskov     12 February 2025
================================================================================================="""

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
coord = ['ATOM      1  N   HIS A   1      49.668  24.248  10.784  1.00 25.00           N\n',
         'ATOM      2  CA  HIS A   1      50.197  24.248  10.436  1.00 16.00           C\n',
         'ATOM      3  C   HIS A   1      49.169  26.701  10.917  1.00 16.00           C\n',
         'ATOM      4  O   HIS A   1      48.241  26.524  11.749  1.00 16.00           O\n',
         'ATOM      5  CB  HIS A   1      51.312  26.524   9.843  1.00 16.00           C\n'
         ]

def sortyzx(record):
    return (record[1],record[2],record[0])

xyz = []
for line in coord:
    line = line.rstrip()
    col = line.split()
    xyz.append( [ col[6], col[7], col[8] ] )	# store coordinates as a list of x, y, z

for atom in xyz:
    print(f'{atom[0]}\t{atom[1]}\t{atom[2]}')

print()
for atom in sorted(xyz):
    print(f'{atom[0]}\t{atom[1]}\t{atom[2]}')

print()
for atom in sorted(xyz, key=sortyzx):
    # atom fields are y, z, x
    print(f'{atom[0]}\t{atom[1]}\t{atom[2]}')

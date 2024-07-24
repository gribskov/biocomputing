"""=================================================================================================
Summarize the output from pilon for each cycle. Sum up each type of error across all contigs and
report as a table.

Input is the STDERR output from the program. Example:
Pilon version 1.24 Thu Jan 28 13:00:45 2021 -0500
Genome: pilon1/vechat_masurca.pilon1.fasta
Fixing snps, indels, gaps, local
Input genome size: 61734324
Scanning BAMs
vechat_masurca.sorted.bam: 24875715 reads, 0 filtered, 21564610 mapped, 20704113 proper, 690328 stray, FR 100% 354+/-72, max 570
Processing contig_233:F:6607:contig_111:F_pilon:1-317126
frags vechat_masurca.sorted.bam: coverage 43
Total Reads: 111913, Coverage: 43, minDepth: 5
Confirmed 316676 of 317126 bases (99.86%)
Corrected 8 snps; 1 ambiguous bases; corrected 3 small insertions totaling 3 bases, 5 small deletions totaling 10 bases
# Attempting to fix local continuity breaks
# fix break: contig_233:F:6607:contig_111:F_pilon:44959-44962 0 -0 +0 NoSolution
# fix break: contig_233:F:6607:contig_111:F_pilon:48517-49451 0 -0 +0 NoSolution
# fix break: contig_233:F:6607:contig_111:F_pilon:66133-66159 0 -0 +0 NoSolution
# fix break: contig_233:F:6607:contig_111:F_pilon:68448-68464 0 -0 +0 NoSolution
# fix break: contig_233:F:6607:contig_111:F_pilon:114855-114898 0 -0 +0 NoSolution
# fix break: contig_233:F:6607:contig_111:F_pilon:143298-143390 0 -0 +0 NoSolution
# fix break: contig_233:F:6607:contig_111:F_pilon:204671-204797 0 -0 +0 NoSolution
# fix break: contig_233:F:6607:contig_111:F_pilon:212287-212339 0 -0 +0 NoSolution TandemRepeat 1066
# fix break: contig_233:F:6607:contig_111:F_pilon:212842-212872 0 -0 +0 NoSolution TandemRepeat 1064
# fix break: contig_233:F:6607:contig_111:F_pilon:268425-268554 0 -0 +0 NoSolution
Fix mismatch: loc=68479 ref=A was=a
Fix mismatch: loc=66427 ref=A was=a
Fix mismatch: loc=65866 ref=A was=a
Fix mismatch: loc=65824 ref=GGGGG was=ggggg
contig_233:F:6607:contig_111:F_pilon:1-317126 log:
Finished processing contig_233:F:6607:contig_111:F_pilon:1-317126
Processing contig_54:F:5278:contig_134:F_pilon:1-123849
frags vechat_masurca.sorted.bam: coverage 44
Total Reads: 43647, Coverage: 44, minDepth: 5
Confirmed 123594 of 123849 bases (99.79%)
Corrected 1 snps; 0 ambiguous bases; corrected 1 small insertions totaling 1 bases, 0 small deletions totaling 0 bases
# Attempting to fix local continuity breaks
# fix break: contig_54:F:5278:contig_134:F_pilon:11897-11980 0 -0 +0 NoSolution
# fix break: contig_54:F:5278:contig_134:F_pilon:23197-23251 0 -0 +0 NoSolution
# fix break: contig_54:F:5278:contig_134:F_pilon:33989-33991 0 -0 +0 NoSolution
# fix break: contig_54:F:5278:contig_134:F_pilon:38426-38495 0 -0 +0 NoSolution TandemRepeat 12
# fix break: contig_54:F:5278:contig_134:F_pilon:42225-42253 0 -0 +0 NoSolution
# fix break: contig_54:F:5278:contig_134:F_pilon:42528-42612 0 -0 +0 NoSolution
# fix break: contig_54:F:5278:contig_134:F_pilon:44484-44533 0 -0 +0 NoSolution
# fix break: contig_54:F:5278:contig_134:F_pilon:44934-44939 0 -0 +0 NoSolution
# fix break: contig_54:F:5278:contig_134:F_pilon:52403-

Michael Gribskov     24 July 2024
================================================================================================="""
import sys


def read_data(filelist):
    """---------------------------------------------------------------------------------------------
    Read in all the data and return a list of hashes, one hash for each cycle. Each cycle hash
    contains
        {'read_n':          read_n,
            'bases_confirmed': bases_confirmed,
            'bases_total':     bases_total,
            'snps':            snps,
            'bases_ambig':     bases_ambig,
            'insert':          insert,
            'bases_insert':    bases_insert,
            'delete':          delete,
            'bases_delete':    bases_delete
        }
        
    :param filelist: list       list of target files from command line split on commas
    :return: dict               a list of summary statistics for each cycle
                                each hash element has a vector of n_cycles values
    ---------------------------------------------------------------------------------------------"""
    cycles = {'read_n':          [],
              'bases_confirmed': [],
              'bases_total':     [],
              'snps':            [],
              'bases_ambig':     [],
              'insert':          [],
              'bases_insert':    [],
              'delete':          [],
              'bases_delete':    []}

    for filename in filelist:
        file = open(filename, 'r')
        cycles['read_n'].append(0)
        cycles['bases_confirmed'].append(0)
        cycles['bases_total'].append(0)
        cycles['snps'].append(0)
        cycles['bases_ambig'].append(0)
        cycles['insert'].append(0)
        cycles['bases_insert'].append(0)
        cycles['delete'].append(0)
        cycles['bases_delete'].append(0)

        for line in file:

            if line.startswith('Total Reads'):
                # total mapped reads to this contig
                # Total Reads: 43647, Coverage: 44, minDepth: 5
                field = line.split(' ')
                cycles['read_n'][-1] += int(field[2])

            elif line.startswith('Confirmed'):
                # confirmed bases
                # Confirmed 123594 of 123849 bases (99.79%)
                field = line.split(' ')
                cycles['bases_confirmed'][-1] += int(field[2])
                cycles['bases_total'][-1] += int(field[4])

            elif line.startswith(''):
                # Corrected 1 snps; 0 ambiguous bases; corrected 1 small insertions totaling 1 bases, 0 small deletions totaling 0 bases
                field = line.split(' ')
                cycles['snps'][-1] += int(field[2])
                cycles['bases_ambig'][-1] += int(field[4])
                cycles['insert'][-1] += int(field[8])
                cycles['bases_insert'][-1] += int(field[12])
                cycles['delete'][-1] += int(field[14])
                cycles['bases_delete'][-1] += int(field[18])

            else:
                # skip all other lines
                continue

        file.close()

    return cycles


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    filelist = sys.argv[1].split(',')
    print(f'{len(filelist)} files to examine')
    cycle = 0
    for filename in filelist:
        print(f'\tcycle {cycle}: {filename}')
        cycle += 1

    cycles = read_data(filelist)

    # when the file is finished save the summary counts
    out = open('pilon_summary,txt', 'w')
    column = ['read_n', 'bases_confirmed', 'bases_total', 'snps', 'bases_ambig',
               'insert', 'bases_insert', 'delete', 'bases_delete']

    legend = {'read_n':          'Number of reads',
              'bases_confirmed': 'Confirmed bases',
              'bases_total':     'Total bases',
              'snps':            'SNPs',
              'bases_ambig':     'Ambiguous bases',
              'insert':          'Insertions',
              'bases_insert':    'Bases in insertions',
              'delete':          'Deletions, ',
              'bases_delete':    'Bases in Deletions'

              }

    c = 0
    for row in cycles:
        stat = column[c]
        # each stat is a row of values for each cycle
        out.write(f'{legend[stat]:22s}' )
        for v in cycles[stat]:
            out.write(f'{v:10d}')

        exit(0)

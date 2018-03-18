"""=================================================================================================
a manual pipeline

Michael Gribskov     March 18  2018
================================================================================================="""
import subprocess as sub

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    project = 'pipe_m1'
    input = ['data/s1_R1.fastq', 's1_R2.fastq']

    # create the project directory
    mkproject = sub.run(['mkdir', project], stdout=sub.PIPE, stderr=sub.PIPE,
                        universal_newlines=True, shell=True)
    prnt('commnd:', mkproject.args())
    print('out|err:', mkproject.stdout, '|', mkproject.stderr)

    exit(0)

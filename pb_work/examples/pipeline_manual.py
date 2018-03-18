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
    input = ['s1_R1.fastq', 's1_R2.fastq']

    # create the project directory
    dir = project
    mkproject = sub.run(['mkdir {}'.format(project)], stdout=sub.PIPE, stderr=sub.PIPE,
                       universal_newlines=True, shell=True)


    # create the subdirectory for the head step
    dir += '/head'
    mkhead = sub.run(['mkdir {}'.format(dir)], stdout=sub.PIPE, stderr=sub.PIPE,
                       universal_newlines=True, shell=True)

    # run the head command
    for filein in input:
        headout = dir + '/' + filein.replace('fastq', 'fastq.head')
        print('head output:', headout)
        head = open( headout, 'w')

        command = 'head -100 {}'.format(filein)
        runhead = sub.call([command ], stdout=head, stderr=sub.PIPE,
                       universal_newlines=True, shell=True)

        head.close()


    exit(0)

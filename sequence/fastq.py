class Fastq:
    """=============================================================================================
    fastq.py

    A simple sequential fastq class
    ============================================================================================="""
    def fastq_read(fh):
        """-----------------------------------------------------------------------------------------
        Read the next 4 lines, a fastq read from the filehandle, fh.
        :param fh: open filehandle
        :return: fq, dict, keys = /title seq sep qual/
        -----------------------------------------------------------------------------------------"""

        fq = {'title': fh.readline().rstrip(),
              'seq': fh.readline().rstrip(),
              'sep': fh.readline().rstrip(),
              'qual': fh.readline().rstrip()}
        # fq['title'] = fh.readline()
        # fq['seq'] = fh.readline()
        # fq['sep'] = fh.readline()
        # fq['qual'] = fh.readline()

        return fq


    def fastq_write(fastq, fh):
        """-----------------------------------------------------------------------------------------

        :param fastq: dictionary with fastq information: keys - title, seq, sep, qual
        :param fh: output filehandle
        :return:
        -----------------------------------------------------------------------------------------"""
        fh.write('{}\n'.format(fastq['title']))
        fh.write('{}\n'.format(fastq['seq']))
        fh.write('{}\n'.format(fastq['sep']))
        fh.write('{}\n'.format(fastq['qual']))

        return


    def sequence_revcomp(seq):
        """----------------------------------------------------------------------------------------
        reverse complement the sequence
        :param seq:
        :return: string, reverse complement of sequence--
        -----------------------------------------------------------------------------------------"""
        return seq.translate(complement)[::-1]


    def quality_rev(qual):
        """-----------------------------------------------------------------------------------------
        Reverse the quality vector to correspond to reverse complement of sequence
        :param qual: list of ascii quality scores
        :return: list (reversed)
        -----------------------------------------------------------------------------------------"""
        return qual[::-1]

# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':


    exit(0)
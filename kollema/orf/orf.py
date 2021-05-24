from sequence.fasta import Fasta


class Orf:
    """=============================================================================================
    Open reading frames. All spans in the reading frame that do not include a stop codon.

    Michael Gribskov     24 May 2021
    ============================================================================================="""

    stop = ('TAG', 'TAA', 'TGA')

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        :param self:
        :return:
        -----------------------------------------------------------------------------------------"""
        self.sequence = ''
        self.rflist = []

    def find(self, direction='+', frame=0):
        """-----------------------------------------------------------------------------------------
        find the open reading frames in a specific frame and direction. For the reverse
        complement, the coordinates are in terms of the reversed sequence

        :param direction: string, '+' or '-'
        :param frame: int, 0 - 2
        :return: int, number of rfs added to self.list
        -----------------------------------------------------------------------------------------"""
        seq = self.sequence
        if direction == '-':
            seq = Fasta.reverseComplement(self.sequence)

        nrf = 0
        pos = frame
        begin = pos

        while pos < len(seq) - 2:
            codon = seq[pos:pos + 3]
            if codon in Orf.stop:
                # end of an ORF
                if pos - begin > 3:
                    nrf += 1
                    self.rflist.append({'direction': direction,
                                        'frame':     frame,
                                        'begin':     begin,
                                        'end':       pos})
                begin = pos + 3

            pos += 3

        if pos - begin > 2:
            nrf += 1
            self.rflist.append({'direction': direction,
                                'frame':     frame,
                                'begin':     begin,
                                'end':       pos})
        return nrf


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    orf = Orf()
    orf.sequence = 'TAAATGATGTGACCCTCACCGTGA'
    print(orf.sequence)

    nrf = 0
    for direction in ('+', '-'):
        s = orf.sequence
        if direction == '-':
            s = Fasta.reverseComplement(orf.sequence)

        for frame in range(3):
            nrf += orf.find(direction=direction, frame=frame)
            print(f'{nrf} reading frames found')

    for i in range(nrf):
        rf = orf.rflist[i]
        s = orf.sequence
        if rf['direction'] == '-':
            s = Fasta.reverseComplement(orf.sequence)
        begin = rf["begin"]
        end = rf["end"]
        print(f'f:{rf["frame"]}{rf["direction"]}\tbegin:{begin:4d}\tend:{end:4d}\t{s[begin:end]}')

    exit(0)

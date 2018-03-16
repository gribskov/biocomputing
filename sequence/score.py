class score():
    """=============================================================================================
    sequence scoring table object
    ============================================================================================="""
    import sys
    def __init__(self, alphabet='ACGT'):
        """-----------------------------------------------------------------------------------------
        score constructor
        -----------------------------------------------------------------------------------------"""
        self.alphabet = ''
        self.a2i = {}
        self.i2a = []
        self.table = []

        alphabetIndex()

    def alphabetIndex(self):
        """-----------------------------------------------------------------------------------------
        Initialze the a2i and i2f for converting between sequence characters and indices
        :return: alphabet size
        -----------------------------------------------------------------------------------------"""
        i = 0
        for char in alphabet:
            if char in self.a2i:
                # duplicate character

            self.a2i[char] = i
            self.i2a.append(char)
            i += 1

        return i

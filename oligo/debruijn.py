####################################################################################################
# for making deBruijn graph examples
####################################################################################################
import sys
import re


class KmerSet():
    """=============================================================================================
    a single kmer
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.text = ''
        k = 0
        self.set = {}

    def add(self, kmer):
        """-----------------------------------------------------------------------------------------
        If unknown add a new kmer to the dict.  ir known, increment count
        :param kmer: string
        :return: int, count of kmer
        -----------------------------------------------------------------------------------------"""
        if kmer in self.set:
            self.set[kmer]['count'] += 1
        else:
            self.set[kmer] = {'count': 1, 'before': [], 'after': []}

        return self.set[kmer]['count']

    def from_text(self, k):
        """-----------------------------------------------------------------------------------------
        populate the kmer set from a string

        :param k: int, length of kmer
        -----------------------------------------------------------------------------------------"""
        text = self.text
        self.k = k

        for i in range(0, len(text) - k + 1):
            kmer = text[i:i + k]
            self.add(kmer)

        return len(self.set.keys())

    def link(self):
        """-----------------------------------------------------------------------------------------
        link each kmer to the previous and following words
        this is simple brute force
        :return:
        -----------------------------------------------------------------------------------------"""
        set = self.set
        kmer = list(set.keys())
        k = self.k

        nlinks = 0
        for i in range(0, len(kmer)):
            k1 = kmer[i]
            suffix = k1[1:]
            prefix = k1[:k - 1]
            for j in range(0, len(kmer)):
                k2 = kmer[j]
                if k2.startswith(suffix):
                    if k1 not in set[k2]['before']:
                        set[k2]['before'].append(k1)
                    if k2 not in set[k1]['after']:
                        set[k1]['after'].append(k2)
                if k2.endswith(prefix):
                    if k1 not in set[k2]['after']:
                        set[k2]['after'].append(k1)
                    if k2 not in set[k1]['before']:
                        set[k1]['before'].append(k2)
                    # if kmer[0].startswith(k2[1:]):
                    #     set[k2]['after'].append(kmer[0])
                nlinks += 1
        # if kmer[0].startswith(kmer[-1][1:]):
        #     set[kmer[-1]]['after'].append(kmer[0])

        return nlinks

    def tips(self):
        """-----------------------------------------------------------------------------------------
        find the left and right tips.  these are the kmers with no befor and after words,
        respectively
        :return:
        -----------------------------------------------------------------------------------------"""
        set = self.set
        left = []
        for kmer in set:
            if len(set[kmer]['before']) == 0:
                left.append(kmer)

        right = []
        for kmer in set:
            if len(set[kmer]['after']) == 0:
                right.append(kmer)

        return left, right


def clean_text(text):
    """---------------------------------------------------------------------------------------------
    remove punction and spaces, force to lower case
    :param text: string
    :return: string, text
    ---------------------------------------------------------------------------------------------"""
    text = text.lower()
    punc = re.compile(r'[^a-z]+')
    text = punc.sub('', text)

    return text


####################################################################################################
# main
####################################################################################################
if __name__ == '__main__':
    text = '''To be, or not to be--that is the question:
        Whether 'tis nobler in the mind to suffer
        The slings and arrows of outrageous fortune
        Or to take arms against a sea of troubles
        And by opposing end them.'''

    text1 = 'To be, or not to be or to be or not'
    text = clean_text(text)
    print(text)

    kmer = KmerSet()
    kmer.text = text
    kmer.from_text(5)
    kmer.link()

    left, right = kmer.tips()
    print('left tips: {}'.format(left))
    print('right tips: {}'.format(right))

    set = kmer.set
    keys = list(set.keys())
    stack = []
    used = {}
    indent = 0
    space = '  '
    if left:
        for kmer in left:
            stack.append([kmer, indent])
    else:
        stack.append([keys[0], indent])

    indent_current = 0
    while stack:
        kmer, indent = stack.pop()

        if indent_current != indent:
            print('\n{}'.format(space * indent), end='')
            indent_current = indent

        if kmer in used:
            print('{}->'.format(kmer), end='')
            indent_current += 1
            continue
        else:
            print('{}'.format(kmer), end='')
            used[kmer] = True

        nafter = len(set[kmer]['after'])
        if nafter == 0:
            # end of chain, outdent
            print('.', end='')
            indent_current -= 1
            continue
        elif nafter == 1:
            # chain continues, don't change indent
            print('-', end='')
        else:
            # fork, indent
            print('-', end='')
            indent = indent_current + 1

        for k in set[kmer]['after']:
            stack.append([k, indent])

    exit(0)

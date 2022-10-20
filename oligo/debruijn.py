####################################################################################################
# for making deBruijn graph examples
####################################################################################################
import sys
import re
import random


class KmerSet():
    """=============================================================================================
    a single kmer
    ============================================================================================="""

    def __init__(self, text=""):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.text = text
        self.reads = []
        self.k = 0
        self.set = {}

    def sample_reads(self, readlen, coverage):
        """-----------------------------------------------------------------------------------------
        sample a list of words of the specified length with at least the specified coverage
        :param readlen: int, length of read to sample
        :param coverage: float, mnimum coverage
        :retur: float, actual coverage
        -----------------------------------------------------------------------------------------"""
        self.reads = []
        text = self.text
        letters = 0
        textlen = len(text)
        while letters / textlen < coverage:
            pos = random.randrange(textlen - readlen)
            self.reads.append(text[pos:pos + readlen])
            letters += readlen

        return letters / textlen

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
        if len(self.reads) == 0:
            self.reads.append(self.text)
        self.k = k

        for text in self.reads:
            for i in range(0, len(text) - k + 1):
                kmer = text[i:i + k]
                self.add(kmer)

        return len(self.set.keys())

    def from_text_random(self, k, n):
        """-----------------------------------------------------------------------------------------
        Randomly sample nword kmer words from the text and
        :return:
        -----------------------------------------------------------------------------------------"""
        if len(self.reads) == 0:
            self.reads.append(self.text)
        self.k = k

        for text in self.reads:
            for _ in range(n):
                pos = random.randrange(len(text) - k)
                kmer = text[pos:pos + k]
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
        find the left and right tips.  these are the kmers with no before and after words,
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

    def list_by_alpha(self):
        """-----------------------------------------------------------------------------------------
        return a string listing the words alphabetically
        :return:
        -----------------------------------------------------------------------------------------"""
        set = self.set
        words = ''
        for kmer in sorted(set):
            words += '{}\t{}'.format(kmer, set[kmer]['count'])
            if 'label' in set[kmer]:
                words += f'\t{set[kmer]["label"]}'
            words += '\n'

        return words

    def clean_text(self):
        """---------------------------------------------------------------------------------------------
        remove punctuation and spaces, force to lower case
        :param text: string
        :return: int, length of text
        ---------------------------------------------------------------------------------------------"""
        text = self.text.lower()
        punc = re.compile(r'[^a-z]+')
        self.text = punc.sub('', text)

        return len(self.text)


####################################################################################################
# main
####################################################################################################
if __name__ == '__main__':
    text = '''To be, or not to be--that is the question:
        Whether 'tis nobler in the mind to suffer
        The slings and arrows of outrageous fortune
        Or to take arms against a sea of troubles
        And by opposing end them.'''

    text1 = 'To be, or not to be--that is the question: Whether tis nobler in the mind to suffer'
    text2 = 'A A T G C G C T A C G T A G G G T A A T A T A A G A C C A'
    champagne_problems = """You booked the night train for a reason
So you could sit there in this hurt
Bustling crowds or silent sleepers
You're not sure which is worse
Because I dropped your hand while dancing
Left you out there standing
Crestfallen on the landing
Champagne problems
Your mom's ring in your pocket
My picture in your wallet
Your heart was glass, I dropped it
Champagne problems
You told your family for a reason
You couldn't keep it in
Your sister splashed out on the bottle
Now no one's celebrating
Dom Perignon, you brought it
No crowd of friends applauded
Your hometown skeptics called it
Champagne problems
You had a speech, you're speechless
Love slipped beyond your reaches
And I couldn't give a reason
Champagne problems"""
# Your Midas touch on the Chevy door
# November flush and your flannel cure
# "This dorm was once a madhouse"
# I made a joke, "Well, it's made for me"
# How evergreen, our group of friends
# Don't think we'll say that word again
# And soon they'll have the nerve to deck the halls
# That we once walked through
# One for the money, two for the show
# I never was ready so I watch you go
# Sometimes you just don't know the answer
# 'Til someone's on their knees and asks you
# "She would've made such a lovely bride
# What a shame she's fucked in the head," they said
# But you'll find the real thing instead
# She'll patch up your tapestry that I shred
# ’tis the damn season
# And hold your hand while dancing
# Never leave you standing
# Crestfallen on the landing
# With champagne problems
# Your mom's ring in your pocket
# Her picture in your wallet
# You won't remember all my
# Champagne problems
# You won't remember all my
# Champagne problems"""


    kmer = KmerSet(text=champagne_problems)
    print("original text\n", kmer.text)
    textlen = kmer.clean_text()
    print("\ncleaned text\n", kmer.text)
    print("\n{} letters".format(textlen))

    actual_coverage = kmer.sample_reads(16, 3)
    print("coverage: {:.2f}".format(actual_coverage))
    for read in sorted(kmer.reads, key=lambda k: kmer.text.index(k)):
        print("\t {}".format(read))

    k = 6
    # kmer.from_text_random(k, int(coverage * textlen / k))
    kmer.from_text(k)
    print("\nkmer words")
    print(kmer.list_by_alpha())
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
    # if left:
    #     for kmer in left:
    #         stack.append([kmer, indent])
    # else:
    #     # no tips, start at the first word
    #     stack.append([keys[0], indent])

    while len(used) < len(set):

        indent = 0
        for start in sorted(set, key=lambda k: (len(set[k]['before']),set[k]['count'])):
            if start in used:
                continue
            stack.append([start,indent])
            break

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
                print('.\n', end='')
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

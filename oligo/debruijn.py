####################################################################################################
# for making deBruijn graph examples
####################################################################################################
import sys
import re
import random
from oligo.kmer import Kmer
# from oligo.debruijn import KmerSet
# import matplotlib.pyplot as plt
import igraph


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
            # if 'label' in set[kmer]:
            #     words += f'\t{set[kmer]["label"]}'
            if set[kmer]["label"] != kmer:
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


def debruijn_condense(kmerset):
    """---------------------------------------------------------------------------------------------
    Condense nodes in the de Bruijn graph that have only one  neighbor upstream and down stream

    :param kmerset: KmerSet object containing the graph
    :return: KmerSet object with the condensed graph
    ---------------------------------------------------------------------------------------------"""
    set = kmerset.set
    ids = list(set.keys())
    k = kmerset.k

    while ids:
        v1 = ids.pop()
        if len(set[v1]['after']) == 1:
            v2 = set[v1]['after'][0]
            if v2 == 'tobe':
                print('tobe')
            if v2 not in set:
                continue
            if len(set[v2]['before']) == 1:

                new = set[v1]['label'] + set[v2]['label'][k - 1:]
                set[v1]['after'] = set[v2]['after']
                set[v1]['label'] = new

                del set[v2]
                if v2 in ids:
                    ids.remove(v2)

    return kmerset


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    text = '''To be, or not to be--that is the question:
        Whether 'tis nobler in the mind to suffer
        The slings and arrows of outrageous fortune
        Or to take arms against a sea of troubles
        And by opposing end them.'''

    # text1 = 'To be, or not to be--that is the question: Whether tis nobler in the mind to suffer'
    text1 = 'To be, or not to be--that is the question'

    text2 = 'A A T G C G C T A C G T A G G G T A A T A T A A G A C C A'

    text3 = '''O Romeo, Romeo, wherefore art thou Romeo?
        Deny thy father and refuse thy name.
        Or if thou wilt not, be but sworn my love
        And I’ll no longer be a Capulet.
        ‘Tis but thy name that is my enemy:
        Thou art thyself, though not a Montague.
        What’s Montague? It is nor hand nor foot
        Nor arm nor face nor any other part
        Belonging to a man. O be some other name.
        What’s in a name? That which we call a rose
        By any other name would smell as sweet;
        So Romeo would, were he not Romeo call’d,
        Retain that dear perfection which he owes
        Without that title. Romeo, doff thy name,
        And for that name, which is no part of thee,
        Take all myself.'''

    text4 = '''Fourscore and seven years ago our fathers brought forth, on this continent, a new nation, conceived in 
    liberty, and dedicated to the proposition that all men are created equal. Now we are engaged in a great civil war, 
    testing whether that nation, or any nation so conceived, and so dedicated, can long endure. We are met on a great 
    battle-field of that war. We have come to dedicate a portion of that field, as a final resting-place for those who 
    here gave their lives, that that nation might live.
    It is altogether fitting and proper that we should do this.
    But, in a larger sense, we cannot dedicate, we cannot consecrate—we cannot hallow—this ground. The brave men, 
    living and dead, who struggled here, have consecrated it far above our poor power to add or detract.
    The world will little note, nor long remember what we say here, but it can never forget what they did here.
    It is for us the living, rather, to be dedicated here to the unfinished work which they who fought here have thus 
    far so nobly advanced. It is rather for us to be here dedicated to the great task remaining before us—that from 
    these honored dead we take increased devotion to that cause for which they here gave the last full measure of 
    devotion—that we here highly resolve that these dead shall not have died in vain—that this nation, under God, shall 
    have a new birth of freedom, and that government of the people, by the people, for the people, shall not perish 
    from the earth.'''

    text5 = '''Fourscore and seven years ago our fathers brought forth, on this continent, a new nation, conceived in 
        liberty, and dedicated to the proposition that all men are created equal. Now we are engaged in a great civil war, 
        testing whether that nation, or any nation so conceived, and so dedicated, can long endure. We are met on a great 
        battle-field of that war. We have come to dedicate a portion of that field, as a final resting-place for those who 
        here gave their lives, that that nation might live.'''

    text6 = '''But soft, what light through yonder window breaks?
    It is the east, and Juliet is the sun.
    Arise, fair sun, and kill the envious moon,
    Who is already sick and pale with grief
    That thou, her maid, art far more fair than she. . . .
    The brightness of her cheek would shame those stars
    As daylight doth a lamp; her eye in heaven
    Would through the airy region stream so bright
    That birds would sing and think it were not night.'''

    text7 = '''ATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTG
TCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGA
TGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGC
AACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTT
GGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAA
CTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGT
GCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTA
ACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACAT
GGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTG
ATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCC
TGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGA
TCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTG
GCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCA
AGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTAC
CGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTG
GTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCA
AGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCA
CAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAG
ATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCAT
CGTAA'''

    black_dog = '''Hey hey mama said the way you move
Gonna make you sweat, gonna make you groove
Ah, ah, child, way you shake that thing
Gonna make you burn, gonna make you sting.
Hey hey baby when you walk that way
Watch your honey drip, can't keep away
Oh yeah, oh yeah, oh, ah, ah
Oh yeah, oh yeah, oh, ah, ah.
I gotta roll, can't stand still
Got a flamin' heart, can't get my fill
Eyes that shine, burnin' red
Dreams of you all through my head
Ah ah, ah ah, ah ah, ah ah, ah ah, ah ah, ahhh
Hey, baby, whoa baby, pretty baby
Darlin' makes 'em do me now
Hey, baby, oh baby, pretty baby
Move me like you're doin' now
Didn't take too long 'fore I found out
What people mean by down and out
Spent my money, took my car
Started tellin' her friends she gonna be a star
I don't know, but I been told
A big-legged woman ain't got no soul
Oh yeah, oh yeah, ah, ah, ah
Oh yeah, oh yeah, ah, ah, ah
All I ask for when I pray
A steady rollin' woman won't come my way
Need a woman gonna hold my hand
Tell me no lies, make me a happy man
Ah ah, ah ah, ah ah, ah ah, ah ah, ah ah, ahhh.
Ah, yeah!'''

    both_sides_now = '''Rows and flows of angel hair
And ice cream castles in the air
And feather canyons everywhere
Looked at clouds that way
But now they only block the sun
They rain and they snow on everyone
So many things I would have done
But clouds got in my way
I've looked at clouds from both sides now
From up and down and still somehow
It's cloud illusions I recall
I really don't know clouds at all'''

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
    Champagne problems"""
    # You told your family for a reason
    # You couldn't keep it in
    # Your sister splashed out on the bottle
    # Now no one's celebrating
    # Dom Perignon, you brought it
    # No crowd of friends applauded
    # Your hometown skeptics called it
    # Champagne problems"""
    # You had a speech, you're speechless
    # Love slipped beyond your reaches
    # And I couldn't give a reason
    # Champagne problems"""
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
    counting_stars = """
    Lately, I've been, I've been losing sleep
Dreaming about the things that we could be
But baby, I've been, I've been praying hard
Said, "No more counting dollars, we'll be counting stars
Yeah, we'll be counting stars"""
# I see this life, like a swinging vine
# Swing my heart across the line
# And in my face is flashing signs
# Seek it out and ye shall find"""

    draggin_my_heart_around="""Baby, you'll come knocking on my front door
Same old line you used to use before
I said, "Hey well, what am I supposed to do?"
I didn't know what I was getting into
Say you've had a little trouble in town
Say you're keeping some demons down
Stop draggin' my
Stop draggin' my 
Stop draggin' my heart around, yeah
It's hard to think about what you wanted
It's hard to think about what you lost
This doesn't have to be the big get even
This doesn't have to be anything at all
I know you really want to tell me goodbye
I know you really want to be your own girl"""
# And baby, you could never look me in the eye
# Oh, when you buckle with the weight of the words
# Stop draggin' my
# Stop draggin' my (stop)
# Stop draggin' my heart around, yeah
# There's people running 'round loose in the world
# Ain't got nothing better to do
# They make a meal of some bright-eyed kid
# You need someone to look after you
# I know you really want to tell me goodbye
# I know you really want to be your own girl
# And baby, you could never look me in the eye
# But why you buckle with the weight of the words (yeah)
# Stop draggin' my
# Stop draggin' my (stop)
# Stop draggin' my heart around, yeah
# Stop draggin' my heart around
# Stop draggin' my heart around (hey, hey, hey)
# Oh, won't you stop draggin' my heart around? (Hey, hey, hey)
# Stop draggin' my heart around"""

    kmer = KmerSet(text=draggin_my_heart_around)
    print("original text\n", kmer.text)
    textlen = kmer.clean_text()
    print("\ncleaned text\n", kmer.text)
    print("\n{} letters".format(textlen))
    kmer.sample_reads(30, 8)
    nread = 0
    for read in kmer.reads:
        nread += 1
        print(f'{nread}\t{read}')
    print()
    kmer.from_text(5)

    kmer.link()
    kset = kmer.set
    for v in kset:
        kset[v]['label'] = v

    print(f'kmers: {len(kmer.list_by_alpha())}')
    print(kmer.list_by_alpha())

    kmer = debruijn_condense(kmer)
    print(kmer.list_by_alpha())
    print(f'condensed kmers: {len(kmer.set)}')
    ncondensed = 0
    lcondensed = 0
    for k in sorted(kmer.set, key=lambda x:len(kmer.set[x]['label']), reverse=True):
        ncondensed += 1
        lcondensed += len(kmer.set[k]["label"])
        print( f'{ncondensed}\t{lcondensed}\t{kmer.set[k]["label"]} ')

    kset = kmer.set

    # make a vertex index

    vidx = {}
    name = []
    label = []
    i = 0
    for v in kset:
        vidx[v] = i
        name.append(v)
        label.append(kset[v]['label'])
        i += 1
    print(label)
    print(f'{i} vertices found')

    # create the edge list
    edges = []
    ecount = 0
    for v in kset:
        for neighbor in kset[v]['after']:
            edges.append((vidx[v], vidx[neighbor]))
            ecount += 1

    print(f'{ecount} edges found')

    # create the igraph graph
    g = igraph.Graph(edges, directed=True)
    g.vs["name"] = name
    g.vs["label"] = label
    # g.vs["label"] = g.vs["name"]
    print(g.get_edgelist())

    # h = g.layout_kamada_kawai(maxiter=50, kkconst=5)
    h = g.layout_fruchterman_reingold(niter=50000, start_temp=1000)
    # h = g.layout_drl()
    # h = g.layout_lgl(maxdelta=75, repulserad=3375000, area=22500, coolexp=1.3, maxiter=5000)
    # h = g.layout_graphopt(niter=10000, node_charge=0.0001, spring_length=10, spring_constant=0.5)
    # h = g.layout_circle()
    # h = g.layout_davidson_harel(maxiter=100, fineiter=50, weight_node_dist=1, weight_border=0, weight_edge_lengths=0.2,
    #                             weight_edge_crossings=2, weight_node_edge_dist=2, cool_fact=0.99)
    style = {'layout':            h,
             'bbox':              (2500, 2500),
             'margin':            600,
             'vertex_label_dist': 5,
             'vertex_label_size': 20,
             'edge_arrow_size':   1,
             'edge_arrow_width':  0.7,
             'edge_width':        2,
             'edge_curved':       True
             }
    igraph.plot(g, "test.pdf", **style)

    exit(0)

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
        If unknown add a new kmer to the dict.  if known, increment count
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

    def from_file(self, filename):
        """-----------------------------------------------------------------------------------------
        open file and read kmers, kmers can be any number per line, space delimited. tokens that are
        only digits are skipped.

        :param filename: string     path to file
        :return: int                number of kmers read
        -----------------------------------------------------------------------------------------"""
        infile = open(filename, 'r')

        lastkmer = ''
        for line in infile:
            token = line.rstrip().split()
            for kmer in token:
                if kmer.isdigit():
                    self.set[lastkmer]['count'] = int(kmer)
                    continue
                self.add(kmer)
                lastkmer = kmer
                l = len(kmer)

        infile.close()
        self.k = l
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
        return a string listing the words alphabetically, if label is defferent than the key, for
        instance, after condensing, show the labell as well

        :return: string     printable list of kmers
        -----------------------------------------------------------------------------------------"""
        set = self.set
        words = ''
        for kmer in sorted(set):
            words += '{}\t{}'.format(kmer, set[kmer]['count'])
            # if 'label' in set[kmer]:
            #     words += f'\t{set[kmer]["label"]}'
            # if set[kmer]["label"] != kmer:
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
        v1 = ids[-1]
        # print(f'v1:{set[v1]}')
        if len(set[v1]['after']) == 1:
            # v1 has only one outgoing edge, may be mergable
            v2 = set[v1]['after'][0]
            # print(f'     v2:{set[v2]}')
            if v2 not in set:
                continue

            if len(set[v2]['before']) == 1:
                # v2 has only one incoming edge, mergable
                newlabel = set[v1]['label'] + set[v2]['label'][k - 1:]
                set[v1]['after'] = set[v2]['after']
                set[v1]['label'] = newlabel

                del set[v2]
                if v2 in ids:
                    ids.remove(v2)
            else:
                # can't extend v1, remove from list
                ids.pop()

        else:
            # if v1 has multiple after nodes, remove from the list
            ids.pop()

    return kmerset


def remove_overlap(kmer):
    """---------------------------------------------------------------------------------------------
    remove the overlap between words by removing the last k-1 characters of each label. Be sure k is
    set in KmerSet object

    :param kmer: KmerSet
    :return:int             sum of lengths of labels
    ---------------------------------------------------------------------------------------------"""
    remove = kmer.k - 1
    total_len = 0
    for k in kmer.set:
        if kmer.set[k]['after']:
            kmer.set[k]['label'] = kmer.set[k]['label'][:-remove]
        total_len += len(kmer.set[k]['label'])

    return total_len

def get_text(name, length=0):
    """---------------------------------------------------------------------------------------------
    A bunch of text stored as a dictionary

    :param name: string     key to look up text
    :param length: int        how many characters to return from the text (after cleaning)
    :return: string         selected text
    ---------------------------------------------------------------------------------------------"""
    library = {
        'hamlet': """
            To be, or not to be--that is the question:
            Whether 'tis nobler in the mind to suffer
            The slings and arrows of outrageous fortune
            Or to take arms against a sea of troubles
            And by opposing end them.""",
        'dna1': 'A A T G C G C T A C G T A G G G T A A T A T A A G A C C A',
        'dna2': """
            ATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTG
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
            CGTAA""",
        'romeo': """
            O Romeo, Romeo, wherefore art thou Romeo?
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
            Take all myself.""",
        'gettysburg': """
            Fourscore and seven years ago our fathers brought forth, on this continent, a new 
            nation, conceived in liberty, and dedicated to the proposition that all men are created 
            equal. Now we are engaged in a great civil war, testing whether that nation, or any 
            nation so conceived, and so dedicated, can long endure. We are met on a great 
            battle-field of that war. We have come to dedicate a portion of that field, as a final 
            resting-place for those who here gave their lives, that that nation might live.
            It is altogether fitting and proper that we should do this. But, in a larger sense, we 
            cannot dedicate, we cannot consecrate — we cannot hallow—this ground. The brave men, 
            living and dead, who struggled here, have consecrated it far above our poor power to 
            add or detract.
            The world will little note, nor long remember what we say here, but it can never forget 
            what they did here. It is for us the living, rather, to be dedicated here to the 
            unfinished work which they who fought here have thus far so nobly advanced. It is 
            rather for us to be here dedicated to the great task remaining before us—that from 
            these honored dead we take increased devotion to that cause for which they here gave 
            the last full measure of devotion—that we here highly resolve that these dead shall 
            not have died in vain—that this nation, under God, shall have a new birth of freedom, 
            and that government of the people, by the people, for the people, shall not perish 
            from the earth.""",
        'butsoft': """
            But soft, what light through yonder window breaks?
            It is the east, and Juliet is the sun.
            Arise, fair sun, and kill the envious moon,
            Who is already sick and pale with grief
            That thou, her maid, art far more fair than she. . . .
            The brightness of her cheek would shame those stars
            As daylight doth a lamp; her eye in heaven
            Would through the airy region stream so bright
            That birds would sing and think it were not night.""",
        'black_dog': """
            Hey hey mama said the way you move
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
            Ah, yeah!""",
        'both_sides_now': """
            Rows and flows of angel hair
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
            I really don't know clouds at all""",
        'champagne_problems': """
            You booked the night train for a reason
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
            You told your family for a reason
            You couldn't keep it in
            Your sister splashed out on the bottle
            Now no one's celebrating
            Dom Perignon, you brought it
            No crowd of friends applauded
            Your hometown skeptics called it
            Champagne problems.
            You had a speech, you're speechless
            Love slipped beyond your reaches
            And I couldn't give a reason
            Your Midas touch on the Chevy door
            November flush and your flannel cure
            "This dorm was once a madhouse"
            I made a joke, "Well, it's made for me"
            How evergreen, our group of friends
            Don't think we'll say that word again
            And soon they'll have the nerve to deck the halls
            That we once walked through
            One for the money, two for the show
            I never was ready so I watch you go
            Sometimes you just don't know the answer
            'Til someone's on their knees and asks you
            "She would've made such a lovely bride
            What a shame she's fucked in the head," they said
            But you'll find the real thing instead
            She'll patch up your tapestry that I shred
            ’tis the damn season
            And hold your hand while dancing
            Never leave you standing
            Crestfallen on the landing
            With champagne problems
            Your mom's ring in your pocket
            Her picture in your wallet
            You won't remember all my
            Champagne problems
            You won't remember all my
            Champagne problems""",
        'counting_stars': """
            Lately, I've been, I've been losing sleep
            Dreaming about the things that we could be
            But baby, I've been, I've been praying hard
            Said, "No more counting dollars, we'll be counting stars
            Yeah, we'll be counting stars
            I see this life, like a swinging vine
            Swing my heart across the line
            And in my face is flashing signs
            Seek it out and ye shall find""",
        'draggin_my_heart_around': """
            Baby, you'll come knocking on my front door
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
            I know you really want to be your own girl
            And baby, you could never look me in the eye
            Oh, when you buckle with the weight of the words
            Stop draggin' my
            Stop draggin' my (stop)
            Stop draggin' my heart around, yeah
            There's people running 'round loose in the world
            Ain't got nothing better to do
            They make a meal of some bright-eyed kid
            You need someone to look after you
            I know you really want to tell me goodbye
            I know you really want to be your own girl
            And baby, you could never look me in the eye
            But why you buckle with the weight of the words (yeah)
            Stop draggin' my
            Stop draggin' my (stop)
            Stop draggin' my heart around, yeah
            Stop draggin' my heart around
            Stop draggin' my heart around (hey, hey, hey)
            Oh, won't you stop draggin' my heart around? (Hey, hey, hey)
            Stop draggin' my heart around""",
        'one_ring': """
            Three Rings for the Elven-kings under the sky,
            Seven for the Dwarf-lords in their halls of stone,
            Nine for Mortal Men doomed to die,
            One for the Dark Lord on his dark throne;
            In the Land of Mordor where the Shadows lie.
            One Ring to rule them all, one Ring to find them,
            One Ring to bring them all, and in the darkness bind them;
            In the Land of Mordor where the Shadows lie.""",
        'gold': """
            All that is gold does not glitter,
            Not all those who wander are lost;
            The old that is strong does not wither,
            Deep roots are not reached by the frost.
    
            From the ashes a fire shall be woken,
            A light from the shadows shall spring;
            Renewed shall be blade that was broken,
            The crownless again shall be king""",
        'only_human': """
            I'm only human, I'm only, I'm only
            I'm only human, human
            Maybe I'm foolish, maybe I'm blind
            Thinking I can see through this and see what's behind
            Got no way to prove it so maybe I'm lying
            But I'm only human after all, I'm only human after all
            Don't put your blame on me, don't put your blame on me
            Take a look in the mirror and what do you see
            Do you see it clearer or are you deceived, in what you believe
            'Cause I'm only human after all, you're only human after all
            Don't put the blame on me
            Don't put your blame on me
            Some people got the real problems
            Some people out of luck
            Some people think I can solve them
            Lord heavens above
            I'm only human after all, I'm only human after all
            Don't put the blame on me
            Don't put the blame on me
            Don't ask my opinion, don't ask me to lie
            Then beg for forgiveness for making you cry, making you cry
            'Cause I'm only human after all, I'm only human after all
            Don't put your blame on me, don't put the blame on me""",
        'declaration': """
        We hold these truths to be self-evident, 
        that all men are created equal, 
        that they are endowed by their Creator with certain unalienable Rights, 
        that among these are Life, Liberty and the pursuit of Happiness.""",
        'holmes': """
        Crime is common. Logic is rare. Therefore it is upon the logic rather than upon the crime that you should dwell"""
        }
    if not length:
        # default length is the entire quote
        length = len(library[name])

    return library[name][:length]


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    generate_kmers = False
    k = 6
    text_len = 0
    text = get_text('holmes')

    if generate_kmers:
        # generate new kmers from text
        kmer = KmerSet(text=text)
        print("original text\n", kmer.text)
        textlen = kmer.clean_text()
        print("\ncleaned text\n", kmer.text)
        print("\n{} letters".format(textlen))
        kmer.sample_reads(30, 4)
        nread = 0
        for read in kmer.reads:
            nread += 1
            print(f'{nread}\t{read}')
        print()
        kmer.from_text(k)
    else:
        # read kmers from list
        kmer = KmerSet()
        kmerfile = 'kmerlist.txt'
        print(f'kmers read from {kmerfile}')
        kmer_n = kmer.from_file(kmerfile)
        print(f'kmers read: {kmer_n}')

    kmer.link()
    kset = kmer.set
    for v in kset:
        kset[v]['label'] = v

    # get total number of kmers
    kmer_total = 0
    for k in kmer.set:
        kmer_total += kmer.set[k]['count']

    print(f'total kmer count: {kmer_total}')
    print(kmer.list_by_alpha())

    kmer = debruijn_condense(kmer)
    print(f'kmers: {len(kmer.list_by_alpha())}')
    print(kmer.list_by_alpha())
    total_len = remove_overlap(kmer)
    print(kmer.list_by_alpha())
    print(f'condensed kmers: {len(kmer.set)} assembly_length: {total_len}')
    half_len = total_len / 2.0
    ncondensed = 0
    lcondensed = 0
    prev_lsum = 0
    prev_l = 0
    n50_notfound = True
    for k in sorted(kmer.set, key=lambda x: len(kmer.set[x]['label']), reverse=True):
        ncondensed += 1
        this_len = len(kmer.set[k]["label"])
        lcondensed += this_len

        if n50_notfound and lcondensed > half_len:
            # found N50
            n50_notfound = False
            frac = (half_len - prev_lsum) / (lcondensed - prev_lsum)
            n50 = prev_l + frac * (this_len - prev_l)
            print(f'\t\t\t--> N50 {half_len:.1f} = {n50:.1f}')

        prev_lsum = lcondensed
        prev_l = this_len

        print(f'{ncondensed}\t{this_len}\t{lcondensed}\t{kmer.set[k]["label"]} ')

    print(f'N50:{n50:.1f}')

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
    print(f'\n{i} vertices found')

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
    # h = g.layout_fruchterman_reingold(niter=50000, start_temp=1000)
    # h = g.layout_drl()
    # h = g.layout_lgl(maxdelta=75, repulserad=3375000, area=22500, coolexp=1.3, maxiter=5000)
    h = g.layout_graphopt(niter=10000, node_charge=0.0001, spring_length=10, spring_constant=0.5)
    # h = g.layout_circle()
    # h = g.layout_davidson_harel(maxiter=100, fineiter=50, weight_node_dist=1, weight_border=0, weight_edge_lengths=0.2,
    #                             weight_edge_crossings=2, weight_node_edge_dist=2, cool_fact=0.99)
    style = {'layout':            h,
             'bbox':              (2500, 2500),
             'margin':            600,
             'vertex_label_dist': 2,
             'vertex_label_size': 20,
             'edge_arrow_size':   1,
             'edge_arrow_width':  0.7,
             'edge_width':        2,
             'edge_curved':       True
             }
    igraph.plot(g, "test.pdf", **style)

    exit(0)

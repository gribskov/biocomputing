"""=================================================================================================


Michael Gribskov     17 April 2022
================================================================================================="""
from kmer import Kmer
from oligo.debruijn import KmerSet
# import matplotlib.pyplot as plt
import igraph


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
    kmer.from_text(6)

    kmer.link()
    kset = kmer.set
    for v in kset:
        kset[v]['label'] = v

    print(kmer.list_by_alpha())
    kmer = debruijn_condense(kmer)
    print(kmer.list_by_alpha())
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
    style = {'layout':     h,
             'bbox':       (2500, 2500),
             'margin':     600,
             'vertex_label_dist': 5,
             'vertex_label_size': 20,
             'edge_arrow_size': 1,
             'edge_arrow_width':0.7,
             'edge_width':2,
             'edge_curved': True
             }
    igraph.plot(g, **style)

    exit(0)

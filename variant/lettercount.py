"""=================================================================================================


Michael Gribskov     08 April 2021
================================================================================================="""

st= 't'
while st:
    count = {}
    st = input('count me?').rstrip().upper()
    st = st.replace(',','.')
    for c in st:

        if c in count:
            count[c.upper()] += 1
        else:
            count[c.upper()] = 1

    for c in count:
        print('{}:{}'.format(c, count[c]))


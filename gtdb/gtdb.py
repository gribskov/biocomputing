"""=================================================================================================
genome taxonomy database data


Michael Gribskov     02 March 2019
================================================================================================="""

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    taxfh = open('1902_bac_metadata_r86.tsv.txt', 'r', encoding="utf-8")

    n = 0
    header = taxfh.readline().rstrip()
    col = header.split('\t')

    taxa = []
    for line in taxfh:
        n += 1
        # print(n)
        field = line.split('\t')
        record = {}
        for i in range(len(field)):
            record[col[i]] = field[i]

        taxa.append(record)
        # if n > 10000:
        #     break

    count = {}
    splittaxa = []
    for record in taxa:
        taxtag = record['gtdb_taxonomy'].split(';')
        if len(taxtag) == 1:
            continue

        bylevel = {}
        for tag in taxtag:
            level, name = tag.split('__')

            bylevel[level] = name

            if not level in count:
                count[level] = {}

            if name in count[level]:
                count[level][name] += 1
            else:
                count[level][name] = 1

        splittaxa.append(bylevel)

    level = 'g'
    minlevel = 60
    selected = []
    for name in count[level]:
        if count[level][name] > minlevel:
            selected.append(name)

    unique = {}
    for record in splittaxa:
        if record[level] in selected:
            ostr = record['p']
            for t in 'cofg':
                ostr += '\t{}'.format(record[t])
            if not ostr in unique:
                unique[ostr] = 1

    nu = 0
    for u in sorted(unique):
        nu += 1
        print('{} {}'.format(nu, u))


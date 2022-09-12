"""=================================================================================================
Read the predicted RNA from the mirDBv6.0 releas and extrct RNAs and targets for a species
================================================================================================="""

if __name__ == '__main__':
    prefix = 'cfa'
    mirdb = open('C:/Users/mgribsko/Dropbox/21dog_dlbcl/miRDB_v6.0_prediction_result.txt', 'r')

    mir_n = 0
    target_n = 0
    unique_n = 0
    target = {}
    unique = {}
    scoremax = 0
    scoremin = 100000
    for line in mirdb:
        # print(target_n, line)
        if not line.startswith(prefix):
            continue

        mir, seq, score = line.rstrip().split()
        score = float(score)

        if score < 80:
            # skip low scores
            continue

        if mir not in target:
            target[mir] = [{'target':seq, 'score':score}]
            mir_n += 1

        else:
            target[mir].append({'target':seq, 'score':score})

        if seq not in unique:
            unique[seq] = 1
            unique_n += 1
        else:
            unique[seq] += 1

        scoremax = max(scoremax, score)
        scoremin = min(scoremin, score)
        target_n += 1

    print(f'{mir_n} miRNA found with prefix {prefix}')
    print(f'{unique_n} unique targets found in {target_n} total targets')
    print(f' Scores in range {scoremin:.2f} to {scoremax:.2f}')

    # histogram of score
    binsize = 1
    hist = [0.0 for _ in range(50)]
    for mir in target:
        for t in target[mir]:
            index = int((t['score'] - 50.0) // binsize)
            hist[index] += 1

    for count in range(len(hist)):
        index = count * binsize + 50.0
        print(f'{index}\t{hist[count]}')

    for seq in unique:
        print(seq)

    exit(0)
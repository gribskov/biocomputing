"""#################################################################################################
ten_pin_bowling.py
$

2025-05-25 gribskov
#################################################################################################"""


def check(fxn, query, result):
    answer = fxn(query)
    print(f'{query}: {answer}, {result}', end=' => ')
    status = 'fail'
    if answer == result:
        status = 'pass'
    print(f'{status}')

    return


def get_balls(score):
    balls = []
    for c in score:
        if c == ' ':
            continue
        elif c == 'X':
            balls.append(10)
        elif c == '/':
            balls.append(10 - balls[-1])
        else:
            balls.append(int(c))

    return balls


def bowling_score(scorestr):
    ball = get_balls(scorestr)
    b = 0
    frame = 0
    score = 0
    nball = 0
    fscore = 0
    for b in range(len(ball)):
        nball += 1
        if frame == 10:
            break
        score += ball[b]
        fscore += ball[b]

        if fscore == 10:
            # strike or spare
            extra = 2 - nball + 1
            for j in range(0, extra):
                score += ball[b + j + 1]
            nball = 2

        if nball == 2:
            nball = 0
            fscore = 0
            frame += 1

    return score


# main

check(bowling_score, '11 11 11 11 11 11 11 11 11 11', 20)
check(bowling_score, 'X X X X X X X X X XXX', 300)
check(bowling_score, '63 44 70 40 9/ 7/ 41 21 62 34', 82)
check(bowling_score, 'X 81 X 05 6/ 9/ X 9/ X 34', 151)
check(bowling_score, '42 51 X 13 35 70 35 X 60 53', 83)

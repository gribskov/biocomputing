"""=================================================================================================
blast.py
Run blast search at ncbi vi blast API

    
================================================================================================="""
import requests
import time

blast = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
program = 'blastp'
database = 'pdb'
query = '''>AAG47671.1 ARV1 [Homo sapiens]
AMGNGGRSGCQYRCIECNQEAKELYRDYNHGVLKITICKSCQKPVDKYIEYDPVIILINAILCKAQAYRHILFNTQINIHGKLYLRWWQLQDSNQNTAPDDLIRYAKEWDF'''

# command = 'Put&PROGRAM={}&DATABASE={}&QUERY={}'.format(program, database, query)
command = {'CMD':         'Put',
           'PROGRAM':     program,
           'DATABASE':    database,
           'QUERY':       query,
           'EMAIL':       'gribskov@purdue.edu',
           }
# print('command:', command)
start = time.asctime( time.localtime(time.time()) )
response = requests.post(blast, command)
print('Blast search started at {}'.format(start))

# get the request ID (RID) and estimated time (RTOE)
infokey = 'QBlastInfoBegin'
infostart = response.text.find(infokey)
ridtag = 'RID = '
ridstart = response.text.find(ridtag, infostart) + len(ridtag)
ridend = response.text.find('\n', ridstart)
rid = response.text[ridstart:ridend]

rtoetag = 'RTOE = '
rtoestart = response.text.find(rtoetag, infostart) + len(rtoetag)
rtoeend = response.text.find('\n', rtoestart)
rtoe = int(response.text[rtoestart:rtoeend])

print('request:{}     estimated time:{}'.format(rid, rtoe))
# rid = '3X4NXWCD016'
# command = 'CMD=Get&&FORMAT_OBJECT=SearchInfo&RID={}'.format(rid)
command = {'CMD':           'Get',
           'FORMAT_OBJECT': 'SearchInfo',
           'RID':           rid
           }

trials_max = 20
interval = max(60,rtoe*15)
result_not_ready = True
trial = 0
print('Polling server for result at {} second intervals'.format(interval))
while result_not_ready:
    # polling loop, don't poll too often, NCBI requests no more than once/minute

    time.sleep(interval)
    current = time.asctime( time.localtime(time.time()) )

    trial += 1
    response = requests.post(blast, command)

    status_key = 'Status='
    status_begin = response.text.find(status_key) + len(status_key)
    status_end = response.text.find('\n', status_begin)
    status = response.text[status_begin:status_end]
    print('\t{}\t{}\t{}'.format(trial, status, current))
    
    if status == 'READY':
        result_not_ready = False
        break

    if trial > trials_max:
        break

if result_not_ready:
    # polling reached maximum number of trials
    print('unable to find result () in {} tries'.format(rid, trials_max))
    exit(0)

# get the final result

print('Blast search completed at {}'.format(current))
command = 'CMD=Get&FORMAT_TYPE=XML&RID={}'.format(rid)
response = requests.get(blast, command)
# print(response.url)
print(response.text)

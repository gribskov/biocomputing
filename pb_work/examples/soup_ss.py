"""=================================================================================================
screen scraping wormbase with beautiful soup

    
================================================================================================="""
# import requests
# from bs4 import BeautifulSoup

# worm = 'http://www.wormbase.org/species/c_elegans/cds/K06C4.1'
#
# response = requests.get(worm)
#
# soup = BeautifulSoup(response.content, 'html.parser')
# # print(soup.prettify())
#
# sequences = soup.findAll('li', {'id': 'sequences'})
# for s in sequences:
#     print(s.prettify())
#     print()

# ="/rest/feed/download/cds/K06C4.1/sequences/K06C4.1
import requests
from bs4 import BeautifulSoup

nws = 'https://forecast.weather.gov/'
wlaf = 'MapClick.php?lat=40.43104000000005&lon=-86.91364999999996'
response = requests.get(nws + wlaf)

soup = BeautifulSoup(response.content, 'html.parser')
current = soup.find('div', {'id': 'current-conditions'})

# <p class="myforecast-current-lrg">32&deg;F</p>
temp = current.find('p', {'class': 'myforecast-current-lrg'})
print('temp:', temp.get_text())
print(temp['class'])


# <div id="current_conditions_detail" class="pull-left">
detail = current.find('div', {'id': 'current_conditions_detail'})

# details table
tr = detail.find_all('tr')
for row in tr:
    td = row.find_all('td')
    print(td[0].get_text(), ':', td[1].get_text().strip())


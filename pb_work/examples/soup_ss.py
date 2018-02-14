"""=================================================================================================
$(PROJECT_NAME}

    
================================================================================================="""
import requests
from bs4 import BeautifulSoup

worm = 'http://www.wormbase.org/species/c_elegans/cds/K06C4.1'

response = requests.get(worm)

soup = BeautifulSoup(response.content, 'html.parser')
print(soup.prettify())

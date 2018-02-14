"""=================================================================================================
beautiful soup

    
================================================================================================="""
import requests
from bs4 import BeautifulSoup

rest = "http://www.wormbase.org/species/c_elegans/cds/K06C4.12"
soup = BeautifulSoup(requests.get(rest).text, 'html.parser')

content = soup.findAll(attrs={"id" : "content"})
for c in content:
    print(c)



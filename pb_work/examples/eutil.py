"""=================================================================================================
eutil.py
retrieve information from NCBI  using eutil REST interface
works 4 March 2024
================================================================================================="""
import requests

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=asthma&usehistory=y
api_key = '645c588aed8ff727e6ed8059d10e7db2ea09'
ncbi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
esearch = 'esearch.fcgi?'

# Get request
query = f'api_key={api_key}&db=protein&term=asthma&usehistory=y'
response = requests.get(ncbi + esearch + query)

print(f'response content:\n{response.content}')
print(f'\nresponse text: {response.text}')

# Post request
params = {'api_key':'645c588aed8ff727e6ed8059d10e7db2ea09',
          'db': 'protein',
          'term': 'asthma',
          'usehistory': 'y'}

esearch = 'esearch.fcgi'
response = requests.post(ncbi + esearch, params)
print(f'\n post response text: {response.text}')

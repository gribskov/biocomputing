"""=================================================================================================
Python program to make  a google map image of a specified location using the Google Static Maps API

Michael Gribskov     06 April 2021
================================================================================================="""
import requests

# my api key (use yours, not mine)
api_key = "AIzaSyBSqLQxLrb-wLoMyqBtj9IIowdRDRIAbuk"
api_url = "https://maps.googleapis.com/maps/api/staticmap?"

# center defines the center of the map, equidistant from all edges of the map.
# zoom defines the zoom level of the map
zoom = 16
center = "40.492,-86.886"
marker = 'color:{}|size:{}|label:{}|{},{}'.format('red', 'mid', 'G', 40.49202,-86.88576)

params = {'center':center, 'zoom':zoom, 'scale':2, 'maptype':'hybrid', 'size':'400x400',
         'markers':marker, 'key':api_key }

r = requests.get(api_url, params)

# open a file to write the map to, 'wb' write binary mode
f = open('map.png', 'wb')
f.write(r.content)
f.close()

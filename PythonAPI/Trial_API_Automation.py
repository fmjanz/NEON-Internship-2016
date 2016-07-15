# -*- coding: utf-8 -*-
"""
Spyder Editor

Trial run of automating API calls. Written in Python 3.5
Frances Janz
National Ecological Observatory Network (NEON)
Created: 5/25/16
Updated: 6/10/16
"""

#Use below for Python 3._
#If using Python 2._ use urllib2
#urlopen is in both libraries
from urllib.request import urlopen
import json
from pandas.io.json import json_normalize

url = "http://api.metagenomics.anl.gov/"

#default is 1; may be updated as new versions become available
version = str('1/')

#Note:NEON data is free to the public and does not require
#an authentication key for the API
#resourcePath = str("annotation/sequence/")#should be a string ending with a ?

#optional string used to filter results
#queryString = str("type=organism&source=RefSeq")

#This should be the MG-RAST ID of the NEON file you
#wish to work with (e.g. mgm4637812.3)
fileID = str("mgm4637823.3?")

#String containing API request for download
downloadRequest = url + version + str("download/") + fileID + str("file=")

callMetaG = url + version + str("metagenome?id=") + fileID

#sample info
getBIOMfile = url + version + str("/profile/") + fileID + str("hit_type=all&type=organism")

#annotatedSequences = url + version + resourcePath + fileID + str("?") + queryString

#make GET request
response = urlopen(downloadRequest).read()
MGdata = json.loads(response.decode('utf8')) #creates JSON dict
json_normalize(MGdata)

response2 = urlopen(callMetaG).read()
MGinfo = json.loads(response.decode('utf8'))

response3 = urlopen(getBIOMfile).read()
BIOMinfo = json.loads(response.decode('utf8'))


import requests
r = requests.get(downloadRequest)
#some dataframe will = r.content

f = open("MGRASTmetaData.txt", "x")
f.write(str(MGdata))
f.close()


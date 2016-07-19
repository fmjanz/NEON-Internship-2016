"""Frances Janz
National Ecological Observatory Netweork (NEON)
Written in Python 3.5 using the Anaconda platform for Windows
Created: Jun 21,2016
Updated: Jun 24,2016

Script for using the MG-RAST API to retrieve function and abundance tables from NEON datasets
See http://api.metagenomics.anl.gov/api.html#matrix for details on the API
Note: NEON data is open access and does not require an API key"""

#Import needed modules. These should all be included with Anaconda
#os is a standard library for setting the working directory
#sys and urllib.error are used here for safe calls and error handling
#urllib.request allows you to work with the urls and actually make the API calls
import urllib.request, urllib.error, json, os, sys
import pandas as pd

#set working directory to wherever you want to store the files
os.chdir("C:\\Users\\fjanz\\Documents\\GitHub\\NEON-Internship-2016\\PythonAPI")

#this is the base url for using the API to make a matrix call
url = "http://api.metagenomics.anl.gov/1/matrix/"

#set to the MG-RAST file number that you want to access
fileID = "4637814.3"

#This retrieves a matrix with functional gene information.
#All of the listed parameters are set to their defaults (see API documentation).
functionCall = url + "function?hit_type=single&group_level=level2&evalue=1&source=Subsystems&result_type=abundance&id=mgm" + fileID 

#This retrieves a matrix with taxonomy info and abundance counts.
#All of the listed parmeters are set to their defaults (see API documentation)
#with the exception of evalue (max e-value cutoff) and source (the reference database used).
organismCall = url + "organism?id=mgm" + fileID + "&evalue=1&source=RefSeq"


#make API call for abundance tables
try:
    opener = urllib.request.urlopen(organismCall)
except urllib.error.HTTPError as e:
    print("Error with HTTP request: %d %s\n%s" % (e.code, e.reason, e.read()))
    sys.exit(255)
    
try:
    opener = urllib.request.urlopen(functionCall)
except urllib.error.HTTPError as e:
    print("Error with HTTP request: %d %s\n%s" % (e.code, e.reason, e.read()))
    sys.exit(255)

matrixData = opener.read()
matrixInfo = json.loads(matrixData.decode('utf8'))

taxonTable = pd.DataFrame.from_dict(matrixInfo,orient='index') #entire matrix file

#You can view taxonTable to get the key names from the JSON object and then
#use whichever one(s) to create the dataset for your file.
#"rows" has the taxon data and "data" the counts
rowData = taxonTable["rows"] #this creates a list

taxonTable2 = pd.DataFrame(data=matrixInfo["data"],index=matrixInfo["rows"])

#normalize takes off leading curly braces, but only works for one key
from pandas.io.json import json_normalize
taxonTable = json_normalize(matrixInfo,['rows'])

#write data to csv file
taxonTable2.to_csv("AnotherTest.csv", mode='x')


#another saving option
with open("TestEverything.txt", 'x') as file_handler:
    for item in rowData:
        file_handler.write("{}\n".format(item))


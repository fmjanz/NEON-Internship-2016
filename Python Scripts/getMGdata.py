# -*- coding: utf-8 -*-
"""
Script for downloading NEON metagenomic data from MG-RAST via the MG-RAST API.
Written in Python 3.5
Frances Janz
National Ecological Observatory Network (NEON)
Created: 6/13/16
Updated: 6/13/16
"""
import os, requests

os.chdir("C:\\Users\\fjanz\\Documents\\PythonAPI")

url = "http://api.metagenomics.anl.gov/1/download/mgm"

#MG-RAST IDs for CPER sites
IDlist = ["4664901.3",
"4637813.3",
"4664858.3",
"4637814.3",
"4637811.3",
"4637817.3",
"4664891.3",
"4637818.3",
"4664866.3",
"4637816.3",
"4637819.3",
"4637823.3",
"4664889.3",
"4637824.3",
"4637822.3",
"4637825.3",
"4637821.3",
"4637828.3",
"4664890.3",
"4637827.3",
"4637829.3",
"4637827.3",
]

for i in range(0,22):
	fileID = IDlist[i]
	downloadRequest = url + fileID + str("?file=700.1") #file=700.1 downloads abundance tables
	r = requests.get(downloadRequest)
	MGdata = r.text #saves data as string
	#write data to new file
	f = open("CPER" + fileID + "genusC.tsv", "x")
	f.write(MGdata)
	f.close()

os.chdir("C:\\Users\\fjanz\\Documents\\MGRAST_DataSets")
	
CompIDList = ["4637813.3", "4637814.3", "4637817.3", "4637818.3", "4637819.3", "4637823.3", "4637824.3", "4637825.3", "4637828.3", "4637829.3", "4637827.3"]

for i in range(0,1):
	
fileID = CompIDList[1]
url2 = "http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeOverview&metagenome=" + fileID + "&action=chart_export&name=organism_genus_hits&file=download." + fileID + ".organism_genus_hits"
r = requests.get(url2)
Cdata = r.text #saves data as string
#f = open(("CPER" + fileID + "pieChart.txt"), "x")
f = open("test1.txt","x")
f.write(Cdata)
f.close()

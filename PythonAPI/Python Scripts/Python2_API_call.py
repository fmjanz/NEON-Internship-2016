#!/usr/bin/env python
'''This script retrieves a metagenome_statistics data structure from the MG-RAST API and
plots a graph using data from the web interface'''

import urllib2, json, sys, os

os.chdir("C:\\Users\\fjanz\\Documents\\GitHub\\NEON-Internship-2016\\metagenome_tables\\Lvl2FunctionalTables\\")

API_URL = "http://api.metagenomics.anl.gov/1"

# assign parameters
metagenomes = ["4637814.3", "4637811.3"]
group_level = "level2"
result_type = "abundance"
source = "Subsystems"

# construct API call 
base_url = API_URL + "/matrix/function"
base_url = base_url + "?group_level=%s&result_type=%s&source=%s&evalue=15&" % (group_level, result_type, source)
base_url = base_url + "&".join( [ "id=mgm%s" % m for m in metagenomes ] ) 

# retrieve the data by sending at HTTP GET request to the MG-RAST API
sys.stderr.write("Retrieving %s\n" % base_url)
try:
    opener = urllib2.urlopen(base_url)
except urllib2.HTTPError, e:
    print "Error with HTTP request: %d %s\n%s" % (e.code, e.reason, e.read())
    sys.exit(255)
opener.addheaders = [('User-agent', 'abundance_matrix.py')]

jsonobject = opener.read()

# convert the data from a JSON structure to a python data type, a dict of dicts.
jsonstructure = json.loads(jsonobject)

# unpack and display the data table
cols = jsonstructure["columns"]
rows = jsonstructure["rows"]
data = jsonstructure["data"]
 
h = { (a, b) : int(c) for (a, b, c) in data } 
sys.stdout.write("Taxon\t") 
for j in range(0, len(cols) ):
    sys.stdout.write(cols[j]["id"] +"\t")
print
for i in range( 0, len(rows)):
    sys.stdout.write(str(rows[i]["id"])+"\t") 
    for j in range( 0, len(cols)):
        try:
            sys.stdout.write(str(h[(i, j)])+"\t" )
        except KeyError:
            sys.stdout.write("0\t")
sys.stdout.write("\n")

h = { (a, b) : int(c) for (a, b, c) in data }
fout = open("functionLvl2Data.tsv",'w')
fout.write("Taxon\t\n") 
for j in range(0, len(cols) ):
    fout.write(cols[j]["id"] +"\t\n")
for i in range( 0, len(rows)):
    fout.write(str(rows[i]["id"])+"\t") 
    for j in range( 0, len(cols)):
        try:
            fout.write(str(h[(i, j)])+"\t"+"\n")
        except KeyError:
            fout.write("0\t")
fout.write("\n")
fout.close()

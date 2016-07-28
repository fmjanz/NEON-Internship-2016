#!/usr/bin/env python
'''This script retrieves a metagenome_statistics data structure from the MG-RAST API and
plots a graph using data from the web interface'''

import urllib, urllib2, json, sys, os
import numpy as np

# convert the data from a JSON structure to a python data type, a dict of dicts.
jsonobject=open(sys.argv[1]).read()
#print "#", sys.argv[1]
jsonstructure = json.loads(jsonobject)

# unpack and display the data table
subsystemsdata = jsonstructure["statistics"]["ontology"]["Subsystems"]
name = jsonstructure["name"]
projectno = jsonstructure["project"][0]
try:
    projectname = jsonstructure["metadata"]["project"]["name"]
except TypeError:
    projectname =""
#ssu = sum(map( int, jsonstructure["statistics"]["source"]["SSU"]["identity"]))
sequence_type = jsonstructure["sequence_type"]
ssufrac = jsonstructure["statistics"]["sequence_stats"]["ratio_reads_rna"]
ssu = jsonstructure["statistics"]["sequence_stats"]["clustered_sequence_count_processed_rna"]
print "#",sys.argv[1] + "\nSSU\t"+ str(ssu)
print "Name\t"+name
print "Projectno\t"+projectno
print "Projectname\t"+projectname
print "sequence_type\t"+sequence_type
print "SSUfrac\t"+ str(ssufrac)
for i in jsonstructure["statistics"]["sequence_stats"].keys():
     print "st."+i + "\t" + jsonstructure["statistics"]["sequence_stats"][i]
for i in jsonstructure["statistics"]["taxonomy"]["domain"]:
     print "tx."+i[0] + "\t" + i[1]
for i in range(len(subsystemsdata)):
    print "ss."+"\t".join(map(str, subsystemsdata[i]))

cols = jsonstructure["statistics"]["qc"]["bp_profile"]["percents"]["columns"]
rows = jsonstructure["statistics"]["qc"]["bp_profile"]["percents"]["data"]
#print "#",sys.argv[1]+"\t"+"\t".join(cols[1:])
for i in range(min(len(rows),50)):
    print "pr."+str(i)+"\t"+str(rows[i][2]+rows[i][3])


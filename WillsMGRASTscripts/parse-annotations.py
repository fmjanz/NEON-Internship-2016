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
data = jsonstructure["statistics"]["ontology"]["Subsystems"]
#ssu = sum(map( int, jsonstructure["statistics"]["source"]["SSU"]["identity"]))
ssufrac = jsonstructure["statistics"]["sequence_stats"]["ratio_reads_rna"]
ssu = jsonstructure["statistics"]["sequence_stats"]["clustered_sequence_count_processed_rna"]
print "#",sys.argv[1] + "\nSSU\t"+ str(ssu)
print "SSUfrac\t"+ str(ssufrac)
for i in range(len(data)):
    print "\t".join(map(str, data[i]))

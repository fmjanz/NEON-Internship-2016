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
cols = jsonstructure["statistics"]["qc"]["bp_profile"]["percents"]["columns"]
rows = jsonstructure["statistics"]["qc"]["bp_profile"]["percents"]["data"]

print "#",sys.argv[1]+"\t"+"\t".join(cols[1:])
for i in range(len(rows)):
    print "\t".join(map(str, rows[i]))

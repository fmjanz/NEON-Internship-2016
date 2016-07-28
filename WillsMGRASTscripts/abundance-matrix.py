#!/usr/bin/env python
'''This script retrieves a metagenome_statistics data structure from the MG-RAST API and
plots a graph using data from the web interface'''

import urllib, sys, os

from mglib import get_auth_token, obj_from_url, biom_to_tab

from optparse import OptionParser

def get_ids(filename):
    mg = []
    for line in open(filename):
        mg.append(line.strip())
    return mg

API_URL = "http://api.metagenomics.anl.gov/1/matrix/"

if __name__ == '__main__':
    usage = "usage: %prog -i <input sequence file> -o <output file>"
    parser = OptionParser(usage)
#    parser.add_option("-i", "--input", dest="input", default=None, help="Input sequence file.")
    parser.add_option("-s", "--source", dest="source", default="RefSeq", help="Annotation source: RefSeq, GenBank, IMG, SEED, TrEMBL, SwissProt, PATRIC, KEG, RDP, Greengenes, LSU, SSU")
    parser.add_option("-g", "--grouplevel", dest="grouplevel", default="domain", help="Grouping level: strain, species, genus, family, order, class, phylum, domain / function, level1, level2, level3")
    parser.add_option("-l", "--list", dest="targetlist", default="", help="Target list (filename).")
#    parser.add_option("-o", "--output", dest="output", default=None, help="Output file.")
    parser.add_option("-i", "--hittype", dest="hittype", default="single", help="Hit type: all, single, lca")
    parser.add_option("-c", "--call", dest="call", default="organism", help="organism or function")
    parser.add_option("-e", "--evalue", dest="evalue", default="1", help="organism or function")
    parser.add_option("-t", "--type", dest="resulttype", default="abundnace", help="Result type: abundnace, evalue, identity, or length")
#    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
    parser.add_option("-k", "--token", dest="token", type="str", help="Auth token")
    parser.add_option("-m", "--metagenomes", dest="metagenomes", default="", type="str", help="Metagenome list")

    (opts, args) = parser.parse_args()
    key = get_auth_token(opts)
# assign parameters
    if not opts.targetlist == "":
        metagenomes = get_ids(opts.targetlist)
    elif not opts.metagenomes == "":
       metagenomes = opts.metagenomes.split(",")
    else:
        metagenomes = ["mgm4447943.3", "mgm4447102.3"]
    group_level = opts.grouplevel
    result_type = "abundance"
    result_call = opts.call
    evalue = opts.evalue
    source = opts.source
    hittype = opts.hittype

# construct API call

    parameters = {"evalue":evalue, "source":source}
    #parameters = {"group_level": group_level, "result_type": result_type, "source":source, "evalue":evalue, "hit_type": hittype}
    base_url = API_URL + result_call + "?" + urllib.urlencode(parameters) +"&" + "&".join(["id=%s" % m for m in metagenomes])

    sys.stderr.write(base_url+"\n")
    jsonstructure = obj_from_url(base_url) # , user_agent="abundance_matrix.py")
    
# unpack and display the data table
    sys.stdout.write("# " + jsonstructure["url"] + "\n")
    cols = jsonstructure["columns"]
    rows = jsonstructure["rows"]
    data = jsonstructure["data"]

    h = {(a, b) : int(c) for (a, b, c) in data}
    sys.stdout.write("MGRID\tMGRID\t")
    for j in range(0, len(cols)):
        sys.stdout.write(cols[j]["id"] +"\t")
    sys.stdout.write("\n")
    sys.stdout.write("MGName\tMGName\t")
    for j in range(0, len(cols)):
        sys.stdout.write(cols[j]["name"] +"\t")
    print
    for i in range(0, len(rows)):
        sys.stdout.write(str(rows[i]["id"])+"\t")
        if result_call == "organism":
            sys.stdout.write(";".join(rows[i]["metadata"]["taxonomy"])+"\t")
        else:
            sys.stdout.write(";".join(rows[i]["metadata"]["ontology"])+"\t")

        for j in range(0, len(cols)):
            try:
                sys.stdout.write(str(h[(i, j)])+"\t")
            except KeyError:
                sys.stdout.write("0\t")
        sys.stdout.write("\n")


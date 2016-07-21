''' Frances Janz
	National Ecological Observatory Network (NEON)
	Created: Jul 19, 2016
	Updated: Jul 19, 2016
	
	This script retrieves functional tables from NEON data stored on MG-RAST.
	Adapted from abundance_table.py from MG-RAST-Tools repository on GitHub.
	
'''
import urllib2, json, sys, os

os.chdir("C:\\Users\\fjanz\\Documents\\GitHub\\NEON-Internship-2016\\metagenome_tables\\Lvl2FunctionalTables\\")

IDlist = ["4637812.3", "4664901.3", "4637813.3", "4664858.3", "4637811.3", "4637814.3", "4637817.3", "4664891.3", "4637818.3", "4664866.3",
"4637816.3", "4637819.3", "4637823.3", "4664889.3", "4637822.3", "4637824.3", "4637821.3", "4637825.3", "4637828.3", "4664890.3", "4637827.3",
"4637829.3", "4637836.3", "4637837.3", "4637835.3", "4637838.3", "4637834.3", "4637839.3", "4637841.3", "4637843.3", "4637842.3", "4637844.3",
"4637840.3", "4637845.3", "4637848.3", "4664912.3", "4664918.3", "4664926.3", "4637849.3", "4664870.3", "4664865.3", "4664925.3", "4664897.3",
"4664910.3", "4664915.3", "4664922.3", "4664893.3", "4664909.3", "4664868.3", "4664892.3", "4664863.3", "4664873.3", "4664920.3", "4664877.3",
"4664886.3", "4664900.3", "4664852.3", "4664914.3", "4637854.3", "4637855.3", "4664884.3", "4664850.3", "4637857.3", "4637858.3", "4637860.3",
"4637861.3", "4637859.3", "4637862.3", "4637864.3", "4637865.3", "4637863.3", "4637866.3", "4637867.3", "4637868.3", "4664894.3", "4664905.3"]

API_URL = "http://api.metagenomics.anl.gov/1"

# assign parameters
metagenomes = IDlist
group_level = "level2"
result_type = "abundance"
source = "Subsystems"
evalue = "2"

# construct API call 
base_url = API_URL + "/matrix/function"

def makeAPIcall (mgID, url):
	# retrieve the data by sending at HTTP GET request to the MG-RAST API
	sys.stderr.write("Retrieving %s\n" % url)
	try:
		opener = urllib2.urlopen(url)
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

	#write data to file
	h = { (a, b) : int(c) for (a, b, c) in data }
	fout = open(mgID + "functionLvl2.tsv",'w')
	fout.write("Taxon\t\n") 
	for j in range(0, len(cols) ):
		fout.write(cols[j]["id"] +"\t\n")
	for k in range( 0, len(rows)):
		fout.write(str(rows[k]["id"])+"\t") 
		for j in range( 0, len(cols)):
			try:
				fout.write(str(h[(k, j)])+"\t"+"\n")
			except KeyError:
				fout.write("0\t")
	fout.write("\n")
	fout.close()


for i in range(0, len(metagenomes)):
	#make call
	full_url = base_url + "?group_level=%s&result_type=%s&source=%s&evalue=%s&" % (group_level, result_type, source, evalue) + "id=mgm" + metagenomes[i]
	makeAPIcall(metagenomes[i],full_url)
	



	

	
	



# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 13:43:35 2016

@author: nrobinson
"""
#Import packages
import mechanize
import cookielib
from bs4 import BeautifulSoup
import html2text
import urllib2
import pandas as pd

#Set login info and urls
username = 'neoninc'
password = 'omicsrus'
login_url = 'http://metagenomics.anl.gov/'
proj_url = 'http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=13948'

#Define functions
#Function to get colum names and data from individual mgrast tables
def extractData (lst,tblName):
    cols = [item[1] for item in lst if item[0] == tblName]
    dat = [item[2] for item in lst if item[0] == tblName]
    return cols, dat 

###############################################################################
#Login
#Set cookies and required browser settings
br = mechanize.Browser()
cj = cookielib.LWPCookieJar()
br.set_cookiejar(cj)

# Browser options
br.set_handle_equiv(True)
br.set_handle_gzip(True)
br.set_handle_redirect(True)
br.set_handle_referer(True)
br.set_handle_robots(False)
br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)

br.addheaders = [('User-agent', 'Chrome')]

# Open main site, where login occurs, and look at forms to make sure that 
# the correct login credentials are applied
br.open(login_url)

#Print forms to take a gander
for f in br.forms():
    print f

br.select_form(nr=0)

# Set user credentials (these are the two forms that were NOT readonly, so user can enter text)
br.form['login'] = username
br.form['password'] = password

# Login
br.submit()

###############################################################################
#Convert html text to object and parse
htmlText=br.open(proj_url).read()

#Extract list of MG-RAST IDs for which to download metadata
#Get IDs (this is a bunch of messy text still)
idList = htmlText.split("id='list_select_preselect_1'")[1].split("input type=")[0]

#Clean up messy text and convert to list of IDs
idList = idList.replace(' value=',''); idList = idList.replace('><','')
idList = idList.replace('~@',' ');idList = idList.replace("'","")
idList = idList.split()

del htmlText    #Cleanup
###############################################################################
#Download metadata files for each ID in idList
#Set main url string for downloading metadata files, with ID to be appended to end for actual download
meta_url = 'http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeOverview&metagenome='

for ind in range(0,len(idList)):
    print 'working on ' + str(idList[ind])
    try:    
        page = meta_url + str(idList[ind])
        #Open page, scrape metadata table for relevant info
        metaText = br.open(page).read()
        #Clean up messy text and convert to list of rows in table
        metaList = metaText.split("id='table_data_0'")[1].split("input type=")[0]
        metaList = metaList.replace(' value=',''); metaList = metaList.replace('>\n<','')
        metaList = metaList.replace('@~',',');metaList = metaList.replace("'","")
        metaList = metaList.replace('NEON, a','NEON a'); metaList = metaList.replace('Boulder, CO, USA','Boulder CO USA')
        metaList = metaList.replace('NEON, Inc.','NEON Inc.'); metaList = metaList.replace('@^',';')
        metaList = metaList.split(','); metaList = [x.split(';') for x in metaList]
        
        #Write to output
        #Get column names and data
        (projCols,projData) = extractData(metaList,'Project')  
        (sampCols,sampData) = extractData(metaList,'Sample')
        (libCols,libData) = extractData(metaList,'Library: metagenome')
        (epCols,epData) = extractData(metaList,'Enviromental Package: soil')
        #Initialize data fromes using first record, then add data from others. Include re-ordering
        # to make sure all data come out in same order in final dataframe
        if ind == 0:
            projDF = pd.DataFrame([index for (elem,index) in sorted(zip(projCols, projData))])
            projDF = projDF.transpose(); projDF.columns = sorted(projCols)
            sampDF = pd.DataFrame([index for (elem,index) in sorted(zip(sampCols, sampData))])
            sampDF = sampDF.transpose(); sampDF.columns = sorted(sampCols)
            libDF = pd.DataFrame([index for (elem,index) in sorted(zip(libCols, libData))])
            libDF = libDF.transpose(); libDF.columns = sorted(libCols)
            epDF = pd.DataFrame([index for (elem,index) in sorted(zip(epCols, epData))])
            epDF = epDF.transpose(); epDF.columns = sorted(epCols)
        else:
            projDF.loc[projDF.shape[0]] = [index for (elem,index) in sorted(zip(projCols, projData))]
            sampDF.loc[sampDF.shape[0]] = [index for (elem,index) in sorted(zip(sampCols, sampData))]
            libDF.loc[libDF.shape[0]] = [index for (elem,index) in sorted(zip(libCols, libData))]
            epDF.loc[epDF.shape[0]] = [index for (elem,index) in sorted(zip(epCols, epData))] 
    except:
        print 'cannot get table for ' + str(idList[ind])



#Reorder columns correctly and write outputs to csvs (later step = combine into one Excel with tabs?)
projDF = projDF[['project_name','mgrast_id','PI_email','PI_firstname','PI_lastname','PI_organization',
'PI_organization_address','PI_organization_country','envo_release', 'project_description']]
projDF.to_csv('C:/Users/nrobinson/Desktop/mgrast_Project.csv',index=False)
sampDF = sampDF[['sample_name','mgrast_id','biome','collection_date','country','env_package','feature','latitude',
'location','longitude', 'material','collection_time','collection_timezone','sample_id']]
sampDF.to_csv('C:/Users/nrobinson/Desktop/mgrast_Samples.csv',index=False)
libDF = libDF[['sample_name','mgrast_id','metagenome_id','metagenome_name','investigation_type','seq_meth',
'file_checksum','file_name','library_name']]
libDF.to_csv('C:/Users/nrobinson/Desktop/mgrast_Library.csv',index=False)
epDF = epDF[['mgrast_id','sample_name','env_package']]
epDF.to_csv('C:/Users/nrobinson/Desktop/mgrast_EP.csv',index=False)

#Question: if there is a missing value somewhere, does my parsing handle it well?
#E.g., if the 1st two elements should be: 
#      'Project;project_name;NEON', 'Project;mgrast_id;mgp13948',
# and there is a missing value, which of the following  comes out? :
#      'Project;project_name;', 'Project;mgrast_id;mgp13948', [THIS IS DESIRED]
#      'Project;project_name;Project', 'mgrast_id;mgp13948;Project', [THIS REPRESENTS A CELL SHIFT- BAD!]


###############################################################################
#Download metadata file for NEON project - this only includes published metagenome and is incomplete
myfile = urllib2.urlopen(meta_url)
output = open('C:/Users/nrobinson/Desktop/forOffline/mgRast/NEON_project_metadata.xlsx', 'wb')
output.write(myfile.read())
output.close()

###############################################################################

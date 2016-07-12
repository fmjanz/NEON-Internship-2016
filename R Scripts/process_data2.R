#Frances Janz
#NEON
#Created: Jun 23, 2016
#Updated: Jul 08, 2016

#Script for processing the metagenome files from the contractor 

#----------------------------------------------------------------------------------------------------------------------
#Load libraries and functions
library('stringr')
library('lessR') #for Sort()
library('vegan') #ecology functions
library("dplyr") 

#merges two files by a specified column
mergeData <- function(fileA, fileB){
  combinedData <- merge(x = fileA, y = fileB, by="Taxonomy", all=TRUE)
  combinedData[is.na(combinedData)] <- 0 #replace NAs with 0s
  return(combinedData)
}

#removes uneeded columns
cleanFile <- function(afile){
  afile <- afile[,-4] #delete phantom empty column
  afile <- afile[,-1] #delete taxon column
  return(afile)
}

#filters out non-target genomic data
cleanTaxon <- function (df){
  df <- filter(df,grepl("virus",df$Taxonomy)==FALSE)
  df <- filter(df,grepl("Eukaryota",df$Taxonomy)==FALSE)
  df <- filter(df,grepl("Archaea",df$Taxonomy)==FALSE)
  df <- filter(df,grepl("unclassified",df$Taxonomy)==FALSE)
  return(df)
}

fill.taxonomy.newMG <- function(x) {
  ### Parses taxon names in abundance table. Specific for 'new' MG-RAST output 
  ### tables with structure "metagenomeID.DB.taxonlevel.csv"
  ### Table should be read in as tab-separated .csv file.
  ### By Lee Stanish June 27, 2016
  
  #tax <- x[,1]
  tax <- x
  alltax <- str_split_fixed(as.character(tax), pattern=";",n=7)
  colnames(alltax) <- c("domain", "phylum", "class", "order", "family","genus", "species")
  alltax <- data.frame(alltax)
  return(alltax)
}

t.otu <- function(x, s, e) {
  ## Written by Lee Stanish
  ## transpose OTU table for ordination analysis##
  ## x=OTU table; s=starting column; e=ending column
  x <- as.data.frame(x) 
  x1 <- t(x[,c(s:e)])      # transpose OTU counts
  colnames(x1) <- x[,1]     #replace OTU ids as column names 
  x1
}

#---------------------------------------------------------------------------------------------------------------------
#Run main script

## Set working directory
if (file.exists(
  "C:/Users/fjanz/Documents")){
  setwd("C:/Users/fjanz/Documents")
}

if (file.exists(
  "/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016")
}


#read in file with ID mapping info
IDfile <- read.csv("SampleID_metadata/ID_mapping_file.csv")

#remove uneeded columns
IDfile <- IDfile[,-1]

if (file.exists(
  "C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")
  }

#tell R where to look for things
if (file.exists(
  "/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata/")
  directory <- "/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata/"
}

files <- list.files(directory, full.names = FALSE)

testFile <- readLines(files[1],2539)

file1 <- read.table(files[1],quote = "", sep='\t',header=TRUE, skip = 2)
file2 <- read.table(files[2],quote = "", sep='\t',header=TRUE, skip = 2)

file1 <- cleanFile(file1)
file2 <- cleanFile(file2)


#combine first two files
MGotu <- as.data.frame(mergeData(file1,file2))

END <- length(files)

#read in and merge remaining files
for (i in 3:END){
  tempfile <- read.table(files[i], quote = "", sep="\t", skip = 2,header = TRUE)
  tempfile <- cleanFile(tempfile)
  MGotu <- as.data.frame(mergeData(MGotu,tempfile))
}

#filter unwanted data
MGotu <- cleanTaxon(MGotu)

#split out taxonomic info
MGtaxa <- fill.taxonomy.newMG(MGotu$Taxonomy)

#transpose columns and rows and organize to prep for diversity analysis
MGotu_flip <- t.otu(x=MGotu,s=2,e=(END + 1))
MGotu_flip <- Sort(by=row.names,data=MGotu_flip)


#repeat process for composite data
if(file.exists(
  "C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")) {
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
}

if(file.exists(
  "/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")) {
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
  directory <- ("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
}

files <- list.files(directory, full.names = FALSE)

file1 <- read.table(files[1],quote = "", sep='\t', skip = 2, header=TRUE)
file2 <- read.table(files[2],quote = "", sep="\t", skip = 2, header = TRUE)
file1 <- cleanFile(file1)
file2 <- cleanFile(file2)
Cotu <- mergeData(file1,file2)

END2 <- length(files)

for (i in 3:END2){
  tempfile <- read.table(files[i], quote = "", sep="\t", skip = 2, header = TRUE)
  tempfile <- cleanFile(tempfile)
  Cotu <- mergeData(Cotu,tempfile)
}

Cotu <- cleanTaxon(Cotu)
Ctaxa <- fill.taxonomy.newMG(Cotu$Taxonomy)
Cotu_flip <- t.otu(x=Cotu,s=2,e=(END2+1))
Cotu_flip <- Sort(by=row.names,data=Cotu_flip)

#----------------------------------------------------------------------------------------------------------------------
##### Diversity analyses on all 38 files for each sample method #####

Cdiv <- (rep(0,END2))
abundanceC <- (rep(0,END2))
for (i in 1:(END2)){
  Cdiv[i] <- diversity(as.numeric(Cotu_flip[i,]))
  abundanceC[i] <- specnumber(Cotu_flip[i,])
}

MGdiv <- (rep(0,END))
abundanceM <- (rep(0,END))
for (i in 1:(END)){
  MGdiv[i] <- diversity(as.numeric(MGotu_flip[i,]))
  abundanceM[i] <- specnumber(MGotu_flip[i,])
}

#convert diversity measures to data frames
Cdiv <- as.data.frame(Cdiv)
MGdiv <- as.data.frame(MGdiv)
abundanceC <- as.data.frame(abundanceC)
abundanceM <- as.data.frame(abundanceM)

#add names for merging
MGplotnames <- rownames(MGotu_flip)
Cplotnames <- rownames(Cotu_flip)
rownames(MGdiv) <- MGplotnames
rownames(Cdiv) <- Cplotnames
rownames(abundanceM) <- MGplotnames
rownames(abundanceC) <- Cplotnames

#-----------------------------------------------------------------------------------------------------------------------
##### combine data ##### 
sampleTypes <- matrix(nrow=END,ncol=2)
sampleTypes[,1] <- rep("Individual",END)
sampleTypes[,2] <- rep("Composite", END2)

#create empty data frame for graphing
divDF <- data.frame(metagenomeID=as.character(),
                    sampleType=as.character(),
                    ShannonIndex=as.numeric(),
                    stringsAsFactors = FALSE)

specRich <- data.frame(metagenomeID=as.character(),
                       sampleType=as.character(),
                       speciesRichness=as.numeric(),
                       stringsAsFactors = FALSE)

e <- length(MGplotnames)*2 #number of samples
j <- 1

#fill graphing data frame
for(i in 1:e){
  if (i <= (e/2)){
    divDF[i,1] <- MGplotnames[i] #create column of sample IDs
    divDF[i,2] <- sampleTypes[i,1] #"Individual"
    divDF[i,3] <- MGdiv[i,1] #diversity measure
    specRich[i,1] <- MGplotnames[i]
    specRich[i,2] <- sampleTypes[i,1]
    specRich[i,3] <- abundanceM[i,1]
  } else {
    divDF[i,1] <- Cplotnames[j]
    divDF[i,2] <- sampleTypes[j,2] #"Composite"
    divDF[i,3] <- Cdiv[j,1]
    specRich[i,1] <- Cplotnames[j]
    specRich[i,2] <- sampleTypes[j,2]
    specRich[i,3] <- abundanceC[j,1]
    j <- j + 1
  }
}

#create columns for easy sorting during graphing
divDF <- merge(divDF,IDfile,by="metagenomeID")
divDF$siteID <- stringr::str_sub(divDF$sampleID,1,4)
divDF$plotID <- stringr::str_sub(divDF$sampleID,1,8)
### Need to create Event_name!!!
ordered <- order(divDF$Event_name)
divDF <- divDF[order(divDF$Event_name),]

specRich <- merge(specRich,IDfile,by="metagenomeID")
specRich$siteID <- stringr::str_sub(specRich$sampleID,1,4)
specRich$plotID <- stringr::str_sub(specRich$sampleID,1,8)
specRich <- specRich[order(specRich$Event_name),]

#merge Individual and Composite OTU tables
merged_OTU <- merge(MGotu,Cotu,by="Taxonomy")

metadata <- data.frame(IDs=divDF$metagenomeID,Events=divDF$Event_name)
#order_metadata <- metadata[order(metadata$Events),]
ordered_OTU <- order(metadata$Events)

#----------------------------------------------------------------------------------------------------------------------
##### Graphs ####


library("ggplot2")
ggplot(divDF,aes(x=Event_name,y=ShannonIndex)) + facet_wrap(~siteID,3,scales="free") +
  geom_bar(aes(fill=sampleType),stat="identity",position="dodge") + scale_fill_manual(values=c("seashell","sandybrown"))

ggplot(specRich,aes(x=plotID,y=speciesRichness)) + facet_wrap(~siteID,3, scales="free") +
  geom_bar(aes(fill=sampleType),stat="identity",position="dodge") + scale_fill_manual(values=c("seashell","sandybrown"))
#----------------------------------------------------------------------------------------------------------------------

##### ANOVAS - Shannon Diversity #####

# model1 <- lm(response ~ predictor_variable, data = data.frame)
# lm(dependent ~ independent)
# lm(y ~ x)
# lm(y ~ x1 + x2 + x3)

##### Models ##### Best Practice: 10 data points per independent/predictor variable
model1 <- lm(ShannonIndex ~ siteID + sampleType, data = divDF)
model2 <- lm(speciesRichness ~ siteID + sampleType, data = specRich)
# model2 <- lm(ShannonIndex ~ siteID + plotID + sampleType, data = divDF)


library("stats")
aShan <- aov(ShannonIndex ~ siteID + sampleType, data = divDF)
aSR <- aov(speciesRichness ~ siteID + sampleType, data=specRich)

##### Summaries
summary(model1)
summary(model2)
summary(aShan)
summary(aSR)

##### line of best fit, x = continuous, y = continuous
x1 <- rnorm(100) # fake predictor 
noise <- rnorm(100, 1,2) # introduce 'natural variation'
y1 <- x1*2.5 + noise # make a response value that is definitely related to x1

####
plot(x = x1, y = y1)
testmod1 <- lm(y1 ~ x1)

summary(testmod1)

##### NMDS ##### 
library(vegan)

?? nms

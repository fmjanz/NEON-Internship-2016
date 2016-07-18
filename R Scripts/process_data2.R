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
library("plyr") #for join()

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
if (file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016")
} else if (file.exists(
  "/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016")
} else {
  print("Error! No such directory.")
}


#read in file with ID mapping info
IDfile <- read.csv("SampleID_metadata/ID_mapping_file.csv")

#remove uneeded columns
IDfile <- IDfile[,-1]

#tell R where to look for things
if (file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")
  } else if (file.exists("/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata/")
  directory <- "/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata/"
  } else {
  print("Error! No such directory.")
}

files <- list.files(directory, full.names = FALSE)

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

#transpose columns and rows and organize to prep for diversity analysis
MGotu_flip <- t.otu(x=MGotu,s=2,e=(END + 1))
MGotu_flip <- Sort(by=row.names,data=MGotu_flip)


#repeat process for composite data
if(file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")) {
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
}else if(file.exists("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")) {
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
  directory <- ("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
} else {
  print("Error! No such directory.")
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

#-----------------------------------------------------------------------------------------------------------------------
##### combine data ##### 

#add names for merging
MGplotnames <- rownames(MGotu_flip)
Cplotnames <- rownames(Cotu_flip)
rownames(MGdiv) <- MGplotnames
rownames(Cdiv) <- Cplotnames
rownames(abundanceM) <- MGplotnames
rownames(abundanceC) <- Cplotnames
sampleTypes <- matrix(nrow=END,ncol=2)
sampleTypes[,1] <- rep("Individual",END)
sampleTypes[,2] <- rep("Composite", END2)

#create empty data frames for graphing
divDF <- data.frame(metagenomeID=as.character(),
                    sampleType=as.character(),
                    ShannonIndex=as.numeric(),
                    stringsAsFactors = FALSE)

specRich <- data.frame(metagenomeID=as.character(),
                       sampleType=as.character(),
                       speciesRichness=as.numeric(),
                       stringsAsFactors = FALSE)

e <- length(MGplotnames)*2 #counter that = total number of samples
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
divDF <- merge(divDF,IDfile,by="metagenomeID") #the fields sampleID and Event_name come with IDfile
divDF$siteID <- stringr::str_sub(divDF$sampleID,1,4)
divDF$plotID <- stringr::str_sub(divDF$sampleID,1,8)

#repeat for species richness object
specRich <- merge(specRich,IDfile,by="metagenomeID")
specRich$siteID <- stringr::str_sub(specRich$sampleID,1,4)
specRich$plotID <- stringr::str_sub(specRich$sampleID,1,8)
specRich <- specRich[order(specRich$Event_name),]

#create and order metadata to use with merged_OTU (see below)
metadata <- divDF[order(divDF$metagenomeID),]
metadata$SpeciesRichness <- specRich$speciesRichness

divDF <- divDF[order(divDF$Event_name),] #sort alpha diversity object

#----------------------------------------------------------------------------------------------------------------------
##### Beta diversity #####
#merge Individual and Composite OTU tables
merged_OTU <- plyr::join(MGotu,Cotu, type="full",match="all")
merged_OTU[is.na(merged_OTU)] <- 0 #replace NAs with 0s

#make taxon names row names
row.names(merged_OTU) <- merged_OTU$Taxonomy
merged_OTU <- merged_OTU[,-1]

#organize beta diversity object
temp_obj <- order(names(merged_OTU))
merged_OTU <- merged_OTU[,temp_obj]

#check that the columns of merged_OTU are in the same order as the ID's in metadata
col <- colnames(merged_OTU)
confirm <- rep(NA,e)
for (i in 1:e){
  if (metadata[i,1] == col[i]){
    confirm[i] <- TRUE
  } else {
    confirm[i] <- FALSE
  }
}

#replaces long sample ID with cleaned version of ID if everything is ordered properly
if (confirm[sample(seq(1,e),1)] == TRUE){
  colnames(merged_OTU) <- metadata$sampleID
}else {
  print("Error! Object names don't match")
}


#save beta diversity object to use with process_OTU_table.r
if (file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016")
}else if (file.exists("/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016")
}

write.csv(merged_OTU,"Metagenome_beta_div_data.csv")
write.csv(metadata,"Metagenome_beta_div_metadata.csv")


#----------------------------------------------------------------------------------------------------------------------
##### Graphs for alpha diversity ####


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

library("stats")
aShan <- aov(ShannonIndex ~ siteID * sampleType, data = divDF) #take out TALL
aSR <- aov(speciesRichness ~ siteID * sampleType, data=specRich)

aShan2 <- aov(ShannonIndex ~ sampleType, data=divDF)
aSR2 <- aov(speciesRichness ~ siteID, data=specRich)

##### Summaries
summary(aShan)
summary(aSR)

TukeyHSD(aShan)
TukeyHSD(aSR2)

#make tables of the overall and the site (?)

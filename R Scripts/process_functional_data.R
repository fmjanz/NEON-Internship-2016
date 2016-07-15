#Frances Janz
#National Ecological Observatory Network
#Adapted from process_data2.R
#Created: Jul 14, 2016
#Updated: Jul 15, 2016

#Script for processing and analyzing functional tables from MG-RAST datasets 


#Load libraries and functions
library('stringr')
library("dplyr")
library('lessR') #for Sort()
library("plyr") #for join()
library('vegan') #diversity()
library("asbio") #evenness()

#removes uneeded columns
cleanFile <- function(afile){
  afile <- afile[,-4] #delete phantom empty column
  afile <- afile[,-1] #delete taxon column
  return(afile)
}

#filters out non-target genomic data
#cleanTaxon <- function (df){
#  df <- filter(df,grepl("virus",df$Taxonomy)==FALSE)
#  df <- filter(df,grepl("Eukaryota",df$Taxonomy)==FALSE)
#  df <- filter(df,grepl("Archaea",df$Taxonomy)==FALSE)
#  df <- filter(df,grepl("unclassified",df$Taxonomy)==FALSE)
#  return(df)
#}

t.otu <- function(x, s, e) {
  ## Written by Lee Stanish
  ## transpose OTU table for ordination analysis##
  ## x=OTU table; s=starting column; e=ending column
  x <- as.data.frame(x) 
  x1 <- t(x[,c(s:e)])      # transpose OTU counts
  colnames(x1) <- x[,1]    #replace OTU ids as column names
  row.names(x1) <- colnames(x[2]) #added line to keep sample ID's attached to datasets -- FJ
  return(x1)
}

#---------------------------------------------------------------------------------------------------------------------
#Run main script

## Set working directory

if(file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/functionTables/MGdata")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/functionTables/MGdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/functionTables/MGdata")
} else if (file.exists("/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata/")
  directory <- "/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata/"
} else {
  print("Error! No such directories.")
}


#read in file with ID mapping info
#IDfile <- read.csv("SampleID_metadata/ID_mapping_file.csv")

#remove uneeded columns
#IDfile <- IDfile[,-1]


files <- list.files(directory, full.names = FALSE)

file1 <- read.table(files[1],quote = "", sep='\t',header=TRUE, skip = 2, fill=TRUE)
file2 <- read.table(files[2],quote = "", sep='\t',header=TRUE, skip = 2, fill=TRUE)

file1 <- cleanFile(file1)
file2 <- cleanFile(file2)

flipped_file1 <- as.data.frame(t.otu(x=file1,s=2,e=2))
flipped_file2 <- as.data.frame(t.otu(x=file2,s=2,e=2))

#combine first two files
MGfunc <- as.data.frame(plyr::join(flipped_file1,flipped_file2,type="full",match="all"))

END <- length(files)

#create vector for sample IDs
sample_names <- rep(NA, END)
sample_names[1] <- row.names(flipped_file1)
sample_names[2] <- row.names(flipped_file2)

#read in and merge remaining files
for (i in 3:END){
  tempfile <- read.table(files[i], quote = "", sep="\t", skip = 2,header = TRUE, fill=TRUE)
  tempfile <- cleanFile(tempfile)
  flip_temp <- as.data.frame(t.otu(x=tempfile,s=2,e=2))
  sample_names[i] <- row.names(flip_temp)
  MGfunc <- as.data.frame(plyr::join(MGfunc,flip_temp,type="full",match="all"))
}

#filter unwanted data
#MGotu <- cleanTaxon(MGotu)

row.names(MGfunc) <- sample_names #attach indentifiers
MGfunc[is.na(MGfunc)] <- 0 #replace NAs with 0s

write.csv(MGfunc,"Joined_Indiv_func_table.csv")

obj <- read.csv("Joined_Indiv_func_table.csv")

#organize to prep for diversity analysis
MGfunc_flip <- Sort(by=row.names,data=MGfunc_flip)


#repeat process for composite data
if(file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/functionTables/Cdata")) {
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/functionTables/Cdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/functionTables/Cdata")
} else if(file.exists("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/functionTables/Cdata")) {
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/functionTables/Cdata")
  directory <- ("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/functionTables/Cdata")
} else {
  print("Error! No such directory.")
}

files <- list.files(directory, full.names = FALSE)

file1 <- read.table(files[1],quote = "", sep='\t', skip = 2, header=TRUE, fill=TRUE)
file2 <- read.table(files[2],quote = "", sep="\t", skip = 2, header = TRUE, fill=TRUE)
file1 <- cleanFile(file1)
file2 <- cleanFile(file2)
Cfunc <- mergeData(file1,file2)

END2 <- length(files)

for (i in 3:END2){
  tempfile <- read.table(files[i], quote = "", sep="\t", skip = 2, header = TRUE, fill=TRUE)
  tempfile <- cleanFile(tempfile)
  Cfunc <- mergeData(Cfunc,tempfile)
}

#Cotu <- cleanTaxon(Cotu)

Cfunc_flip <- t.otu(x=Cotu,s=2,e=(END2+1))
Cfunc_flip <- Sort(by=row.names,data=Cotu_flip)

#----------------------------------------------------------------------------------------------------------------------
##### Diversity analyses on all 38 files for each sample method #####

Cdiv <- (rep(0,END2))
evennessC <- (rep(0,END2))
for (i in 1:(END2)){
  Cdiv[i] <- diversity(as.numeric(Cfunc_flip[i,]))
  evennessC[i] <- evenness(Cfunc_flip[i,]) #change this!!!
}

MGdiv <- (rep(0,END))
evennessM <- (rep(0,END))
for (i in 1:(END)){
  MGdiv[i] <- diversity(as.numeric(MGfunc_flip[i,]))
  evennessM[i] <- evenness(MGfunc_flip[i,]) #change this!!!
}

#convert diversity measures to data frames
Cdiv <- as.data.frame(Cdiv)
MGdiv <- as.data.frame(MGdiv)
evennessC <- as.data.frame(evennessC)
evennessM <- as.data.frame(evennessM)

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
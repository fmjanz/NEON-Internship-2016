#Frances Janz
#National Ecological Observatory Network (NEON)
#Created: Jul 21, 2016
#Updated : Jul 21, 2016

#This script processes functional tables from MG-RAST for diversity analyses
#Adapted from process_data2.R


library('stringr')
library('lessR') #for Sort()
library("plyr") #for join()
library('vegan') #ecology functions

#standardize metagenome IDs for mapping; takes a vector; uses "stringr"
appendMGM <- function (IDs){
  e <- length(IDs)
  for (i in 1:e){
    IDs[i] <- str_c("mgm",as.character(IDs[i]))
  }
  return(IDs)
}

#removes uneeded columns and rows
cleanFile <- function(afile){
  afile <- afile[,-3] #delete phantom empty column
  afile <- afile[-1,] #delete header
  
  #reassign metagenome ID to name of 2nd column
  mgID <- afile[1,1]
  afile <- afile[-1,]
  colnames(afile) <- c("Function",mgID)
  return(afile)
}

#merges two files by a specified column
mergeData <- function(fileA, fileB){
  combinedData <- merge(x = fileA, y = fileB, by="Function", all=TRUE)
  combinedData[is.na(combinedData)] <- 0 #replace NAs with 0s
  return(combinedData)
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

#trades MG-RAST IDs for event names
replaceIDs <- function (IDframe,dataset){
  rnames <- row.names(dataset)
  events <- rep(0,length(dataset[,1]))
  j <- 1
  
  #fills an object with the sample IDs
  for (i in 1:length(IDframe[,1])){
    if (IDframe$mgrastID[i] == rnames[j] && j <= 38){
      events[j] <- IDframe$metagenomeID[i]
      j <- j + 1
    }else {
      next
    }
  }
  
  row.names(dataset) <- events
  return(dataset)
}

#---------------------------------------------------------------------------------------------------------------------
##### Run main script #####

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
IDfile <- read.csv("SampleID_metadata/ID_mapping_file.csv", stringsAsFactors = FALSE)

#standardize metagenome IDs
IDfile[,1] <- appendMGM(IDfile[,1])

#tell R where to look for things
if (file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/MGdata")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/MGdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/MGdata")
} else if (file.exists("/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/MGdata/")
  directory <- "/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/MGdata/"
} else {
  print("Error! No such directory.")
}

files <- list.files(directory, full.names = FALSE)

file1 <- read.csv(files[1],quote = "", sep='\t',header=FALSE, stringsAsFactors = FALSE, fill = TRUE)
file2 <- read.csv(files[2],quote = "", sep='\t',header=FALSE, stringsAsFactors = FALSE, fill = TRUE)

file1 <- cleanFile(file1)
file2 <- cleanFile(file2)


#combine first two files
MGdata <- as.data.frame(mergeData(file1,file2))

END <- length(files)

#read in and merge remaining files
for (i in 3:END){
  tempfile <- read.csv(files[i], quote = "", sep="\t",header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  tempfile <- cleanFile(tempfile)
  MGdata <- as.data.frame(mergeData(MGdata,tempfile))
}

#transpose columns and rows and organize to prep for diversity analysis
MGdata_flip <- t.otu(x=MGdata,s=2,e=(END + 1))
MGdata_flip <- Sort(by=row.names,data=MGdata_flip) #sort IDs in ascending order for next 
MGdata_flip <- replaceIDs(IDfile,MGdata_flip) #replace metagenome IDs with samples events

##### repeat process for composite data ####

if (file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/Cdata")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/Cdata")
  directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/Cdata")
} else if (file.exists("/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/Cdata/")
  directory <- "/Users/lstanish/Github/NEON-Internship-2016/metagenome_tables/Lvl2FunctionalTables/Cdata/"
} else {
  print("Error! No such directory.")
}

files <- list.files(directory, full.names = FALSE)
file1 <- read.csv(files[1],quote = "", sep='\t',header=FALSE, stringsAsFactors = FALSE, fill = TRUE)
file2 <- read.csv(files[2],quote = "", sep='\t',header=FALSE, stringsAsFactors = FALSE, fill = TRUE)
file1 <- cleanFile(file1)
file2 <- cleanFile(file2)

Cdata <- as.data.frame(mergeData(file1,file2))

END2 <- length(files)

#read in and merge remaining files
for (i in 3:END2){
  tempfile <- read.csv(files[i], quote = "", sep="\t",header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  tempfile <- cleanFile(tempfile)
  Cdata <- as.data.frame(mergeData(Cdata,tempfile))
}

Cdata_flip <- t.otu(x=Cdata,s=2,e=(END2 + 1))
Cdata_flip <- Sort(by=row.names,data=Cdata_flip) 
Cdata_flip <- replaceIDs(IDfile,Cdata_flip)

#----------------------------------------------------------------------------------------------------------------------
##### Diversity analyses on all 38 files for each sample method #####

Cdiv <- (rep(0,END2))
evennessC <- (rep(0,END2))

for (i in 1:(END2)){
  Cdiv[i] <- diversity(as.numeric(Cdata_flip[i,]))
  abundance <- specnumber(Cdata_flip[i,])
  evennessC[i] <- Cdiv[i]/log(abundance)
}

MGdiv <- (rep(0,END))
evennessM <- (rep(0,END))


for (i in 1:(END)){
  MGdiv[i] <- diversity(as.numeric(MGdata_flip[i,]))
  abundance <- specnumber(MGdata_flip[i,])
  evennessM[i] <-  MGdiv[i]/log(abundance)
}

#convert diversity measures to data frames
Cdiv <- as.data.frame(Cdiv)
MGdiv <- as.data.frame(MGdiv)
evennessC <- as.data.frame(evennessC)
evennessM <- as.data.frame(evennessM)

#-----------------------------------------------------------------------------------------------------------------------
##### combine data ##### 

#add names for merging
MGplotnames <- rownames(MGdata_flip)
Cplotnames <- rownames(Cdata_flip)
rownames(MGdiv) <- MGplotnames
rownames(Cdiv) <- Cplotnames
rownames(evennessM) <- MGplotnames
rownames(evennessC) <- Cplotnames
sampleTypes <- matrix(nrow=END,ncol=2)
sampleTypes[,1] <- rep("Individual",END)
sampleTypes[,2] <- rep("Composite", END2)

#create empty data frames for graphing
divDF <- data.frame(metagenomeID=as.character(),
                    sampleType=as.character(),
                    ShannonIndex=as.numeric(),
                    stringsAsFactors = FALSE)

even <- data.frame(metagenomeID=as.character(),
                       sampleType=as.character(),
                       Evenness=as.numeric(),
                       stringsAsFactors = FALSE)

e <- length(MGplotnames)*2 #counter that = total number of samples
j <- 1

#fill graphing data frame
for(i in 1:e){
  if (i <= (e/2)){
    divDF[i,1] <- MGplotnames[i] #create column of sample IDs
    divDF[i,2] <- sampleTypes[i,1] #"Individual"
    divDF[i,3] <- MGdiv[i,1] #diversity measure
    even[i,1] <- MGplotnames[i]
    even[i,2] <- sampleTypes[i,1]
    even[i,3] <- evennessM[i,1]
  } else {
    divDF[i,1] <- Cplotnames[j]
    divDF[i,2] <- sampleTypes[j,2] #"Composite"
    divDF[i,3] <- Cdiv[j,1]
    even[i,1] <- Cplotnames[j]
    even[i,2] <- sampleTypes[j,2]
    even[i,3] <- evennessC[j,1]
    j <- j + 1
  }
}

#create columns for easy sorting during graphing
divDF <- merge(divDF,IDfile,by="metagenomeID") #the fields sampleID and Event_name come with IDfile
divDF$siteID <- stringr::str_sub(divDF$sampleID,1,4)
divDF$plotID <- stringr::str_sub(divDF$sampleID,1,8)

#repeat for species richness object
even <- merge(even,IDfile,by="metagenomeID")
even$siteID <- stringr::str_sub(even$sampleID,1,4)
even$plotID <- stringr::str_sub(even$sampleID,1,8)
even <- even[order(even$Event_name),]

#create and order metadata to use with merged_table (see below)
metadata <- divDF[order(divDF$mgrastID),]

#metadata object must be created before this step or the ordering will be off!
divDF <- divDF[order(divDF$Event_name),] #sort alpha diversity object

#----------------------------------------------------------------------------------------------------------------------
##### Beta diversity #####
#merge Individual and Composite OTU tables
merged_table <- plyr::join(MGdata,Cdata, type="full",match="all")
merged_table[is.na(merged_table)] <- 0 #replace NAs with 0s

#make taxon names row names
row.names(merged_table) <- merged_table$Function
merged_table <- merged_table[,-1]

#organize beta diversity object
temp_obj <- order(names(merged_table))
merged_table <- merged_table[,temp_obj]

#check that the columns of merged_table are in the same order as the ID's in metadata
col <- colnames(merged_table)
confirm <- rep(NA,e)
for (i in 1:e){
  if (metadata[i,4] == col[i]){
    confirm[i] <- TRUE
  } else {
    confirm[i] <- FALSE
  }
}

#replaces long sample ID with cleaned version of ID if everything is ordered properly
if (confirm[sample(seq(1,e),1)] == TRUE){
  colnames(merged_table) <- metadata$sampleID
}else {
  print("Error! Object names don't match")
}


#save beta diversity object to use with process_OTU_table.r
if (file.exists("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016")){
  setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016")
}else if (file.exists("/Users/lstanish/Github/NEON-Internship-2016")){
  setwd("/Users/lstanish/Github/NEON-Internship-2016")
}

write.csv(merged_table,"Functional_beta_div_data.csv")
write.csv(metadata,"Functional_beta_div_metadata.csv")


#----------------------------------------------------------------------------------------------------------------------
##### Graphs for alpha diversity ####

library("ggplot2")

plotIDs <- stringr::str_sub(row.names(MGdiv),1,8)

#plot Shannon
ggplot(divDF,aes(x=Event_name,y=ShannonIndex)) + facet_wrap(~siteID,3,scales="free") +
  geom_bar(aes(fill=sampleType),stat="identity",position="dodge") +
  scale_fill_manual(values=c("lightsteelblue2","dodgerblue4")) +
  theme(plot.title=element_text(face="bold"), axis.text.x= element_text(angle=45, hjust = 0.75, vjust = 0.75, size=5.5)) +
  xlab(NULL) + ylab(NULL) + ylim(0,6.5)

#plot Evenness
ggplot(even,aes(x=plotID,y=Evenness)) + facet_wrap(~siteID,3, scales="free") +
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
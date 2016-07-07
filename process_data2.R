#Frances Janz
#NEON
#Created: Jun 23, 2016
#Updated: Jul 05, 2016

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

#filters out non-microbial genomic data
cleanTaxon <- function (df){
  df <- filter(df,grepl("virus",df$Taxonomy)==FALSE)
  df <- filter(df,grepl("Eukaryota",df$Taxonomy)==FALSE)
  df <- filter(df,grepl("Eukaryota",df$Taxonomy)==FALSE)
  df <- filter(df,grepl("unclassified",df$Taxonomy)==FALSE)
  return(df)
}

fill.taxonomy.newMG <- function(x) {
  ### Parses taxon names in abundance table. Specific for 'new' MG-RAST output 
  ### tables with structure "metagenomeID.DB.taxonlevel.csv"
  ### Table should be read in as tab-separated .csv file.
  ### By Lee Stanish June 27, 2016
  
  tax <- x[,1]
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

#creates 6 side-by-side graphs for data; specify colors as strings
makeGraphs <- function(color1,color2,color3,color4){
  leg <- c("Individual", "Composite") #legend names
  par(mfrow=c(2,3))
  sites <- c("CPER","DSNY","HARV","OSBS","STER","TALL")
  
  for (i in 1:(length(sites)-1)){
    barplot(divDF$ShannonIndex[metadataDF$siteIDs==sites[i]],width=0.45,beside = FALSE,ylim=c(5.65,6.2),xlab=sites[i],ylab="Shannon Index",
            horiz=F,xpd=F, las=1, col=c("lavender","seashell"),cex.lab = 1.5, cex.main = 1.4)
  }
  
  barplot(divDF$ShannonIndex[metadataDF$siteIDs==sites[6]],width=0.45,beside = TRUE,ylim=c(5.65,6.2),xlab=sites[6],ylab="Shannon Index",
          legend.text = leg, horiz=F,xpd=F, las=1, col=c(color1,color2),cex.lab = 1.5, cex.main = 1.4)
  
  for (i in 1:(length(sites)-1)){
    barplot(abundance[metadataDF$siteIDs==sites[i]],width=0.45,beside = TRUE,ylim=c(1590,1760),xlab=sites[i],ylab="Species Richness",
            horiz=F,xpd=F, las=1, col=c(color3,color4),cex.lab = 1.5, cex.main = 1.4)
  }
  
  barplot(abundance[metadataDF$siteIDs==sites[6]],width=0.45,beside = TRUE,ylim=c(1590,1760),xlab=sites[6],ylab="Species Richness",
          legend.text = leg, horiz=F,xpd=F, las=1, col=c(color3,color4),cex.lab = 1.5, cex.main = 1.4)
}

#---------------------------------------------------------------------------------------------------------------------
#Run main script
setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")

#tell R where to look for things
directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/MGdata")
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
MGtaxa <- fill.taxonomy.newMG(MGotu)

#transpose columns and rows and organize to prep for diversity analysis
MGotu_flip <- t.otu(x=MGotu,s=2,e=END)
MGotu_flip <- Sort(by=row.names,data=MGotu_flip)


#repeat process for composite data
setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
directory <- ("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/metagenome_tables/OTUtables/Cdata")
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
Ctaxa <- fill.taxonomy.newMG(Cotu)
Cotu_flip <- t.otu(x=Cotu,s=2,e=END2)
Cotu_flip <- Sort(by=row.names,data=Cotu_flip)

#----------------------------------------------------------------------------------------------------------------------
#Diversity analyses on all 40 files for each sample method

Cdiv <- (rep(0,END2-1))
abundanceC <- (rep(0,END2-1))
for (i in 1:(END2-1)){
  Cdiv[i] <- diversity(as.numeric(Cotu_flip[i,]))
  abundanceC[i] <- specnumber(Cotu_flip[i,])
}

MGdiv <- (rep(0,END-1))
abundanceM <- (rep(0,END-1))
for (i in 1:(END-1)){
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
#combine data
div <- as.matrix(cbind(MGdiv,Cdiv))
abundance <- as.matrix(cbind(abundanceM,abundanceC))
plotnames <- rownames(div)

sampleTypes <- matrix(nrow=37,ncol=2)
sampleTypes[,1] <- rep("Individual",37)
sampleTypes[,2] <- rep("Composite", 37)

#create empty data frame for graphing
divDF <- data.frame(SiteID=as.character(),
                    SampleType=as.character(),
                    ShannonIndex=as.numeric(),
                    stringsAsFactors = FALSE)

e <- length(plotnames)*2 #number of samples
j <- 1

#fill graphing data frame
for(i in 1:e){
  if (i <= (e/2)){
    divDF[i,1] <- MGplotnames[i]
    divDF[i,2] <- sampleTypes[i,1]
    divDF[i,3] <- MGdiv[i,1]
  } else {
    divDF[i,1] <- Cplotnames[j]
    divDF[i,2] <- sampleTypes[j,2]
    divDF[i,3] <- Cdiv[j,1]
    j <- j + 1
  }
}

#----------------------------------------------------------------------------------------------------------------------
#Graphs

#create metadata object for sorting graphs by site ID

siteIDs <- rep(NA, length(plotnames))
for (i in 1:length(plotnames)){
  siteIDs[i] <- substr(plotnames[i],1,4)
}
metadata <- cbind(plotnames,siteIDs)
metadataDF <- as.data.frame(metadata)


makeGraphs("wheat3","whitesmoke","sandybrown","seashell")



getMeans <- function (divNum){
  siteMeans <- rep(0,6)
  siteMeans[1] <- mean(divNum[1:11])
  siteMeans[2] <- mean(divNum[12:23])
  siteMeans[3] <- mean(divNum[24:26])
  siteMeans[4] <- mean(divNum[27:32])
  siteMeans[5] <- mean(divNum[33:38])
  siteMeans[6] <- mean(divNum[39])
  return(siteMeans)
}

MGmean <- getMeans(MGdiv)
Cmean <- getMeans(Cdiv)
allMeans <- rbind(MGmean,Cmean)

barplot(allMeans,width=0.5,beside = TRUE,xlim=c(0,12),ylim=c(5.7,7.0),ylab="Shannon Index",
        legend.text = leg, horiz=F,xpd=F, las=1, col=c("maroon","gray"),cex.lab = 1.5, cex.main = 1.4)
axis(side=1,labels= F,at=seq(0,6))


xnum <- seq(1,84)
plot(allMeans,main="Shannon",xlab = "Sites",
     ylab="Shannon Index",type='b')
par(new=T)
plot(Cdiv[metadataDF$siteIDs=="TALL"], col="orange",axes = F,type='b',pch=2)

MGdf <- as.data.frame(MGotu_flip)

xnum <- seq(1,11)
ynum <- seq(5,7)
qplot(x=div[metadataDF$siteIDs=="CPER"])

barplot(divDF$ShannonIndex[metadataDF$siteIDs=="CPER"],width=0.45,beside = TRUE,ylim=c(5.6,6.25),ylab="Shannon Index",
        legend.text = leg, horiz=F,xpd=F, las=1, col=c("tomato","whitesmoke"),cex.lab = 1.5, cex.main = 1.4)

library("ggplot2")
ggplot(divDF,aes(x=SiteID,y=ShannonIndex, fill=SampleType,color=SampleType)) +
  geom_bar(stat="identity",position="dodge")

#fix IDs!!!
for (i in 1:e){
  for (j in 1:7){
    temp <- unlist(str_split(divDF[,1],'[.]', n=8))
    split_IDs[i,j] <- temp[j]
  }
}



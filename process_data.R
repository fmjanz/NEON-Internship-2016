#Frances Janz
#NEON
#Created: 6/23/16
#Updated: 6/29/16

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
  alltax <- str_split_fixed(as.character(tax), pattern=";", 
                            n=7)
  colnames(alltax) <- c("domain", "phylum", "class", "order", "family",
                        "genus", "species")
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

#read in and merge remaining files
for (i in 3:length(files)){
  tempfile <- read.table(files[i], quote = "", sep="\t", skip = 2,header = TRUE)
  tempfile <- cleanFile(tempfile)
  MGotu <- as.data.frame(mergeData(MGotu,tempfile))
}

MGotu <- cleanTaxon(MGotu)

#split out taxonomic info
MGtaxa <- fill.taxonomy.newMG(MGotu)

#transpose columns and rows to prep for diversity analysis
MGotu_flip <- t.otu(x=MGotu,s=2,e=39)

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

for (i in 3:length(files)){
  tempfile <- read.table(files[i], quote = "", sep="\t", skip = 2, header = TRUE)
  tempfile <- cleanFile(tempfile)
  Cotu <- mergeData(Cotu,tempfile)
}

Cotu <- cleanTaxon(Cotu)
Ctaxa <- fill.taxonomy.newMG(Cotu)
Cotu_flip <- t.otu(x=Cotu,s=2,e=39)
Cotu_flip <- Sort(by=row.names,data=Cotu_flip)

#----------------------------------------------------------------------------------------------------------------------
#Diversity analyses on all 40 files for each sample method

Cdiv <- rep(0,38)
abundanceC <- rep(0,38)
for (i in 1:38){
  Cdiv[i] <- diversity(as.numeric(Cotu_flip[i,]))
  abundanceC[i] <- specnumber(Cotu_flip[i,])
}

MGdiv <- rep(0,38)
abundanceM <- rep(0,38)
for (i in 1:38){
  MGdiv[i] <- diversity(as.numeric(MGotu_flip[i,]))
  abundanceM[i] <- specnumber(MGotu_flip[i,])
}

div <- rbind(MGdiv,Cdiv) #combine data

#----------------------------------------------------------------------------------------------------------------------
#Graphs

plotnames <- rownames(MGotu_flip)
siteIDs <- rep(NA, 39)

for (i in 1:length(plotnames)){
  siteIDs[i] <- substr(plotnames[i],1,4)
}

metadata <- cbind(plotnames,siteIDs)

metadataDF <- as.data.frame(metadata)
#CPERids <- match("CPER",siteIDs)

colnames(div) <- plotnames

leg <- c("Individual", "Composite") #legend names

 #fix margins

par(mfrow=c(2,3))
#par(mar=c(3,3,1,1))
barplot(div[metadataDF$siteIDs=="CPER"],width=0.45,beside = TRUE,ylim=c(5.6,6.4),ylab="Shannon Index",
                legend.text = leg, horiz=F,xpd=F, las=1, col=c("blue","cyan"),cex.lab = 1.5, cex.main = 1.4)
axis(side=1,labels= F,at=seq(1,11)) #prep for labels
barplot(div[metadataDF$siteIDs=="DSNY"],width=0.45,beside = TRUE,ylim=c(5.6,6.4),ylab="Shannon Index",
        legend.text = leg, horiz=F,xpd=F, las=1, col=c("purple","gold"),cex.lab = 1.5, cex.main = 1.4)
axis(side=1,labels= F,at=seq(1,12))
barplot(div[metadataDF$siteIDs=="HARV"],width=0.75,beside = TRUE,ylim=c(5.6,6.4),ylab="Shannon Index",
         legend.text=leg,horiz=F,xpd=F, las=1, col=c("maroon","tan"),cex.lab = 1, cex.main = 1)
axis(side=1,labels= F,at=seq(1,3))
barplot(div[metadataDF$siteIDs=="OSBS"],width=0.45,beside = TRUE,ylim=c(5.6,6.4),ylab="Shannon Index",
        legend.text = leg, horiz=F,xpd=F, las=1, col=c("green","dark gray"),cex.lab = 1.5, cex.main = 1.4)
axis(side=1,labels= F,at=seq(1,6))
barplot(div[metadataDF$siteIDs=="STER"],width=0.45,beside = TRUE,ylim=c(5.6,6.4),ylab="Shannon Index",
        legend.text = leg, horiz=F,xpd=F, las=1, col=c("magenta","maroon"),cex.lab = 1.5, cex.main = 1.4)
axis(side=1,labels= F,at=seq(1,6))
barplot(div[metadataDF$siteIDs=="TALL"],width=0.15,beside = TRUE,ylim=c(5.6,6.4),ylab="Shannon Index",
        legend.text = leg, horiz=F,xpd=F, las=1, col=c("yellow","dark green"),cex.lab = 1.5, cex.main = 1.4)
axis(side=1,labels= F,at=seq(0,1))

#text(seq(0, 22, by=2), par("usr")[3] - 0.2, srt = 90, pos = 1, xpd = TRUE) #use Sample ID as x labels

barplot(div[metadataDF$siteIDs=="HARV"],beside = TRUE,col=c("maroon","tan"),ylim=c(5.0,6.35),xpd=F,
        legend.text=leg)


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






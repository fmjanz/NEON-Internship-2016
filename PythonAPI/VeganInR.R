#Frances Janz
#National Ecological Observatory Network (NEON)
#Created: 6/13/16
#Updated: 6/13/16

#Script for finding species richness

install.packages("vegan") #ecology functions
install.packages("reshape") #re-formats matrices
install.packages("plyr")
library("vegan")
library("reshape")
library("plyr")


setwd("C:/Users/fjanz/Documents/MGRAST_DataSets")

directory <- "C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016/MGRAST_DataSets"

files <- list.files(directory, full.names = TRUE)
 file1 <- read.table("CPER4367823.3genusC.tsv",sep='',fill=TRUE,quote="")

#read in and merge user specified tsv files
for (i in 1:length(files)){
  
  file2 <- read.table(files[i], sep='\t', fill = TRUE, quote="")
  
  x <- file2[,1]
  
  temp2 <- data.frame(file2[,2],row.names = x, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
  
  flippedfile <- t(temp2) #switch columns and rows
  
  #name <- paste("new",as.character(i),".tsv",sep="")
  
  #write.table(flippedfile,name,sep="\t",na='0')
  
  file1 <- rbind.fill(as.data.frame(flippedfile),as.data.frame(file1)) #merge data
}


for (i in 2:3525){
  div[i] <- diversity(file1[,i])
}

data(BCI)
#H <- diversity(BCI)

#str(BCI)

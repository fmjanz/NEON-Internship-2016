#Frances Janz
#National Ecological Observatory Network
#Created: 6/1/2016
#Updated: 6/2/2016

#Script for filtering C out of sample names for data comparison of metagenomic data

install.packages("stringr")
install.packages("dplyr")
library("stringr")
library("dplyr")

setwd("C:/Users/fjanz/Documents/Merge_CSV_Files")

file <- read.csv("Merged_Marker_and_MG_Data.csv")

file$DataType <- NA

#collect sample ID's from csv
initial_IDs <- as.character(file$sample_name)
sample_IDs <- as.character(file$sampleID)

END <- length(initial_IDs)

for (i in 1:END){
  if (!is.na(file[i,2])){
    file[i,29] <- "MG"
  } else{
    file[i,29] <- '16S'
  }
}

f <- rep(NA,END)

split_IDs <- data.frame(f,row.names=NULL, check.rows = FALSE,
                        check.names = FALSE, fix.empty.names = FALSE,
                        stringsAsFactors = FALSE)

#splits IDs at separators ('-')
for (i in 1:END){
  for (j in 1:7){
    temp <- unlist(str_split(initial_IDs[i],'[-]', n=7))
    split_IDs[i,j] <- temp[j]
  }
}

#add identifying columns
for (i in 1:END){
  split_IDs[i,7] <- sample_IDs[i] #16S sample IDs
  split_IDs[i,8] <- file[i,29] #MG or 16S
}

#Identifies composite samples and removes C from sample name
for (i in 1:END){
  if (split_IDs[i,4] == 'C'){
    split_IDs[i,6] <- "Composite"
    split_IDs[i,4] <- 'Y'
    split_IDs[i,5] <- split_IDs[i,3]
    split_IDs[i,3] <- 'X'
  } else
    split_IDs[i,6] <- NA
}

MG <- data.frame(Col1=character(),
                 Col2=character(),
                 Col3=character(),
                 Col4=character(),
                 Col5=character(),
                 Col6=character(),
                 Col7=character(),
                 Col8=character(),
                 stringsAsFactors = FALSE)

marker <- data.frame(Col1=character(),
                     Col2=character(),
                     Col3=character(),
                     Col4=character(),
                     Col5=character(),
                     Col6=character(),
                     Col7=character(),
                     Col8=character(),
                     stringsAsFactors = FALSE)
u <- 1
v <- 1

#split marker and MG data
for (i in 1:END){
  if (split_IDs[i,8] == "MG"){
    MG[u,] <- split_IDs[i,]
    u <- u + 1
  } else{
    marker[v,] <- split_IDs[i,]
    v <- v + 1
  }
}

#remove uneeded columns
MG <- MG[,-7]
MG <- MG[,-4]
MG <- MG[,-3]
marker <- marker[,-6]
marker <- marker[,-4]
marker <- marker[,-3]

colnames(MG) <- c('location','horizon',"date","composite","DataType")
colnames(marker) <- c('location','horizon',"date","sampleID","DataType")

matches <- data.frame(Col1=character(),
                      Col2=character(),
                      Col3=character(),
                      Col4=character(),
                      Col5=character(),
                      Col6=character(),
                      Col7=character(),
                      Col8=character(),
                      stringsAsFactors = FALSE)

m <- length(MG$location)
n <- length(marker$location)
j <- 1
k <- 1

#sorts out matching records (very slow!)
for (i in 1:m){
  for (j in 1:n){
    if ((MG[i,1] == marker[j,1]) && (MG[i,2] == marker[j,2]) && (MG[i,3] == marker[j,3])){
      matches[k,] <- marker[j,]
      k <- k + 1
    }else {
      next
    }
  }
  
  matches[k,] <- MG[i,]
  k <- k + 1
}

#remove phantom duplicated columns
matches <- matches[,-8]
matches <- matches[,-7]
matches <- matches[,-6]

e <- length(matches$Col1)

#recombine sample IDs
for (i in 1:e){
  matches[i,1] <- str_c(matches[i,1],matches[i,2],matches[i,3],sep='-')
}

#removes duplicate columns
matches <- matches[,-3]
matches <- matches[,-2]

colnames(matches) <- c('sample_name','sampleID_16S',"dataType")

write.csv(matches,file="Matched_MG_and_16S_Data.csv",row.names = FALSE)
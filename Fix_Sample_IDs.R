#Frances Janz
#National Ecological Observatory Network (NEON)
#Created: 5/26/16
#Updated: 5/31/16

#Script to standardize soil sample ID's across files

install.packages("plyr") #dataframe manipulation
library("plyr")
install.packages("stringr") #string manipulation
library("stringr")

setwd("C:/Users/fjanz/Documents/Merge_CSV_Files")

file1 <- read.csv("Final_NEON_16S_MappingFile_r.csv")
file2 <- read.csv("neon_16s_mapping_file_r.csv")
file3 <- read.csv("MGRAST_metadata_sample.csv")
file3 <- file3[,-1] #remove uncessary rownames

correctProblems <- function(problemFrame){
  n <- nrow(problemFrame)
  
  newFrame <- data.frame(Col1=character(),
                           Col2=character(),
                           Col3=character(),
                           Col4=character(),
                           Col5=character(),
                           Col6=character(),
                           stringsAsFactors = FALSE)
  
  for (i in 1:n){
    if ((problemFrame[i,1] != '') && (problemFrame[i,6] == 'M' || problemFrame[i,6] == 'O')){
      newFrame[i,1] <- str_c(problemFrame[i,1],'_',problemFrame[i,2]) #site Id
      newFrame[i,2] <- problemFrame[i,6] #horizon
      newFrame[i,3] <- problemFrame[i,4] #x coordinates
      newFrame[i,4] <- problemFrame[i,5] #y coordinates
      newFrame[i,5] <- problemFrame[i,7] #date
      newFrame[i,6] <- problemFrame[i,8] #old ID
    }else if ((problemFrame[i,1] != '') && (problemFrame[i,7] == 'M' || problemFrame[i,7] == 'O') && (!is.na(problemFrame[i,7]))){
      newFrame[i,1] <- str_c(problemFrame[i,1],'_',problemFrame[i,2]) #site Id
      newFrame[i,2] <- problemFrame[i,7] #horizon
      newFrame[i,3] <- problemFrame[i,4] #x coordinates
      newFrame[i,4] <- problemFrame[i,5] #y coordinates
      newFrame[i,5] <- problemFrame[i,6] #date
      newFrame[i,6] <- problemFrame[i,8] #old ID
    }else if ((problemFrame[i,1] == '') && (problemFrame[i,7] == 'M' || problemFrame[i,7] == 'O')){
      newFrame[i,1] <- str_c(problemFrame[i,2],'_',problemFrame[i,3]) #site ID
      newFrame[i,2] <- problemFrame[i,7] #horizon
      newFrame[i,3] <- problemFrame[i,5] #x coodinates
      newFrame[i,4] <- problemFrame[i,6] #y coordinates
      newFrame[i,6] <- problemFrame[i,8] #old ID
      newFrame[i,5] <- "DateMissing"
    } else if ((problemFrame[i,1] != '') && (problemFrame[i,4] == 'M' || problemFrame[i,4] == 'O')){
      newFrame[i,1] <- str_c(problemFrame[i,1],'_',problemFrame[i,2]) #site Id
      newFrame[i,2] <- problemFrame[i,4] #horizon
      newFrame[i,3] <- problemFrame[i,5] #x coordinates
      newFrame[i,4] <- problemFrame[i,6] #y coordinates
      newFrame[i,5] <- "DateMissing"
      newFrame[i,6] <- problemFrame[i,8] #old ID
    } else {
      newFrame[i,1] <- "Error"
      newFrame[i,2] <- "Error"
      newFrame[i,3] <- "Error"
      newFrame[i,4] <- "Error"
      newFrame[i,5] <- "Error"
      newFrame[i,6] <- problemFrame[i,8] #old ID
    }
  }
  return(newFrame)
}

#creates a sample_name column with sample IDs in correct format
standardizeSamples <- function (file){
  #collect sample ID's from csv
  initial_IDs <- as.character(file$sampleID)
  
  
  END <- length(initial_IDs)
  f <- rep(NA,END)
  split_IDs <- data.frame(f,row.names=NULL, check.rows = FALSE,
                          check.names = FALSE, fix.empty.names = FALSE,
                          stringsAsFactors = FALSE)
  
  #splits IDs at separators ('.')
  for (i in 1:END){
    for (j in 1:7){
      temp <- unlist(str_split(initial_IDs[i],'[.]', n=8))
      split_IDs[i,j] <- temp[j]
    }
  }
  
  split_IDs$oldIDs <- initial_IDs
  fixedIDs <- data.frame (f,row.names=NULL, check.rows = FALSE,
                          check.names = FALSE, fix.empty.names = FALSE,
                          stringsAsFactors = FALSE)
  
  #create empty data frame
  problemIDs <- data.frame(Col1=character(),
                           Col2=character(),
                           Col3=character(),
                           Col4=character(),
                           Col5=character(),
                           Col6=character(),
                           Col7=character(),
                           Col8=character(),
                           stringsAsFactors = FALSE)
  
  u <- 1 #placeholder value for following for loop
  v <- 1
  
  #Filters on first two info pieces (location and horizon)
  for (i in 1:END) {
    if ((split_IDs[i,3] == "M" || split_IDs[i,3] == "O") && (is.na(split_IDs[i,7]))){
      fixedIDs[v,1] <- str_c(split_IDs[i,1],'_',split_IDs[i,2],sep='')
      fixedIDs[v,2] <- split_IDs[i,3]
      fixedIDs[v,3] <- split_IDs[i,4]
      fixedIDs[v,4] <- split_IDs[i,5]
      fixedIDs[v,5] <- split_IDs[i,6]
      fixedIDs[v,6] <- split_IDs[i,8]
      v <- v + 1
    }else if((split_IDs[i,4] == "M" || split_IDs[i,4] == "O") && (!is.na(split_IDs[i,7]))){
      fixedIDs[v,1] <- str_c(split_IDs[i,1],'_',split_IDs[i,2],sep='')
      fixedIDs[v,2] <- split_IDs[i,4]
      fixedIDs[v,3] <- split_IDs[i,5]
      fixedIDs[v,4] <- split_IDs[i,6]
      fixedIDs[v,5] <- split_IDs[i,7]
      fixedIDs[v,6] <- split_IDs[i,8]
      v <- v + 1
    }else { #sorts out IDs with info in the wrong place
      problemIDs[u,] <- split_IDs[i,]
      u <- u + 1
    }
  }
  
  fixedProblems <- correctProblems(problemIDs)
  colnames(fixedProblems) <- c('V1','V2','V3','V4','V5','V6') #make col names the same for rbind
  colnames(fixedIDs) <- c('V1','V2','V3','V4','V5','V6')
  fixedIDs <- rbind(fixedProblems,fixedIDs)
  
  e <- head(which(is.na(fixedIDs[1])),1) - 1 #find end of data frame containing IDs
  d <- rep(NA,e)
  finalIDs <- data.frame (d,row.names=NULL, check.rows = FALSE,
                          check.names = FALSE, fix.empty.names = FALSE,
                          stringsAsFactors = FALSE)
  
  #put everything back together with separators ("-")
  for (i in 1:e){
    finalIDs[i,1] <- str_c(fixedIDs[i,1],fixedIDs[i,2],fixedIDs[i,3],fixedIDs[i,4], fixedIDs[i,5],sep='-')
    finalIDs[i,2] <- as.character(fixedIDs[i,6])
  }
  colnames(finalIDs) <- c("newID","sampleID")
  
  return(finalIDs)
}

firstMarker <- standardizeSamples(file1)
secondMarker <- standardizeSamples(file2)

#creates new column in original file to store corrected ID's. If corrected ID does not exist, NA is placed in cell
insertNewIDs <- function (dataframe, XFiles){
  temp <- XFiles$sampleID
  e <- length(temp)
  for (i in 1:e){
    j <- match(XFiles$sampleID[i],dataframe$sampleID)
    if (isTRUE(identical(as.character(XFiles$sampleID[i]), dataframe$sampleID[j]))){
      XFiles$sample_name[i] <- dataframe$newID[j]
    }else {
      XFiles$sample_name[i] <- NA  
    }
  }
  return(XFiles)
}

newfile1 <- insertNewIDs(firstMarker,file1)
newfile2 <- insertNewIDs(secondMarker,file2)

#merge files
mergedMarkerFiles <- join(newfile1,newfile2,type="full",match="all")
mergedFiles <- rbind.fill(file3,mergedMarkerFiles)

write.csv(mergedFiles,file="Merged_Marker_and_MG_Data.csv",row.names=FALSE,na='')

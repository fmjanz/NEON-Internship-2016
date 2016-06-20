##### SETUP #####
#install.packages(stringr)
#install.packages(plyr)

library(stringr)
library(plyr)

# point to your folder with all the files
directory <- "C:/Users/cflagg/Documents/GitHub/NEON-Internship-2016/taxonTables"

## OPTION 1: if all of your files are in one folder, just set the working directory
## R only automatically looks at files in a few folders all the time, including your 'working directory'
setwd(directory)

# list all of the files in your folder
files <- list.files(directory, full.names = FALSE)

# check if it works -- open first file
head(read.table(files[1], header = TRUE, fill = TRUE, quote = "", sep="\t"))

## OPTION 2: if your files are in multiple folders, point to the 'root' folder and grab the full file path
# files <- list.files(directory, full.names = TRUE)

###### CUSTOM FUNCTIONS ##### 
# the file list here should point to "fileGrab1" for this specific script
multiCombine <- function(input, ply = ldply){
  ply(input, function(x){
    t <- read.table(x, header=TRUE, sep="\t",fill = TRUE, quote = "") # read the tsv
    ## see what x is
    print(x)
    ## append part of the fileName as a column ## FIX THIS ##
    t$siteID <- stringr::str_sub(x, 1,4)
    ## add mgrast ID from file name
    t$mgrastID <- stringr::str_sub(x,1,10)
    colnames(t) <- c("species", "count", "siteID", "mgrastID")
    t1 <- rbind(t) # rbind it to a temporary variable
    return(t1) # return the full variable
  }
  )
}



# OPTION 2: a re-usable function - outputType can be 'vector' or 'list' 
## use if you need files in more than one folder
#fileParser <- function(directoryList, outputType = "vector"){
  # need to dive into each folder, decide if it's a .csv or .xlsx, then read the rows
  fullFilenameList <- list() # initialize a list object to populate
  counter = 0
  for (folder in directoryList){
    # browser() # this is the interactive debugger, its scope is global & local
    # list all files in a particular folder
    folder_content <- list.files(folder, full.names = TRUE)
    for (fileN in folder_content){
      # iterate the file count for tracking and for indexing the list
      counter = counter + 1
      fullFilenameList[[counter]] <- fileN
      print(fileN)
    }
  }
  # what kind of output should be returned?
  ifelse(outputType == "vector", 
         # return a vector
         return(unlist(fullFilenameList)),
         # return a list
         return(fullFilenameList))
}

##### EXECUTE: Option 1 ##### 
## if using setwd()
storeIt <- multiCombine(files)


##### EXECUTE: Option 2 ##### 
## parse the file paths into a vector
#fileList = fileParser(fileList)

## grep files that you want, if there is more than one type of file in the folder
#fileList_trim <- grep(pattern=".csv", x = fileList)

# merge the individual files, as specified by fileList_trim, into one large data frame
#mergedFiles <- multiCombine(fileList_trim)


##### Subset New Dataset

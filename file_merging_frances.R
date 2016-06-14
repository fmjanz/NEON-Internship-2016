##### SETUP #####
install.packages(stringr)
install.packages(plyr)

library(stringr)
library(plyr)

# point to your folder with all the files
directory <- "C:/francesFiles"

# list all of the files in your folder
files <- list.files(directory, full.names = TRUE)

###### CUSTOM FUNCTIONS ##### 
# the file list here should point to "fileGrab1" for this specific script
multiCombine <- function(input, ply = ldply){
  ply(input, function(x){
    t <- read.csv(x, header=TRUE, sep=",",stringsAsFactors = FALSE) # read the csv
    ## append the fileName as a column ## FIX THIS ##
    ## str_sub
    t1 <- rbind(t) # rbind it to a temporary variable
    return(t1) # return the full variable
  }
  )
}

# a re-usable function - outputType can be 'vector' or 'list'
fileParser <- function(directoryList, outputType = "vector"){
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

##### EXECUTE ##### 
## parse the file paths into a vector
fileList = fileParser(fileList)

## grep files that you want, if there is more than one type of file in the folder
fileList_trim <- grep(pattern=".csv", x = fileList)

# merge the individual files, as specified by fileList_trim, into one large data frame
mergedFiles <- multiCombine(fileList_trim)
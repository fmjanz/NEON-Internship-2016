## Code for compiling a list of files in a target directory for automating script
## running.

#local file locations
setwd("yourdirectory")
file.list <- list.files(path=".", pattern="yourpattern")

for(i in 1:length(file.list)) {
  temp <- read.csv(file.list[i],header = TRUE, sep=",", stringsAsFactors = FALSE) # read in csv
# insert your code here  
}
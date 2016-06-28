setwd("/Users/lstanish/Documents/Data_sequences/Metagenomics/MGRAST_annotation_tables/Run201506_set2/")
x.old <- read.csv("mgm4637810.3.refseq.level6.csv", sep='\t')
x.new <- read.csv("mgm4637809.3.RS.species.csv", sep='\t', skip=2)


fill.taxonomy.oldMG <- function(x) {
  ### Parses taxon names in abundance table. Specific for 'old' MG-RAST output 
  ### tables with structure "metagenomeID.database.level.csv"
  ### Table should be read in as tab-separated .csv file.
  ### By Lee Stanish June 27, 2016
  library(stringr)
  
  tax.p <- str_split_fixed(as.character(x$RefSeq.species), pattern=";", 
                n=8)
  alltax <- tax.p[,2:8]
  colnames(alltax) <- c("domain", "phylum", "class", "order", "family",
                    "genus", "species")
  alltax <- data.frame(alltax)
  return(alltax)
  }
  

fill.taxonomy.newMG <- function(x) {
  ### Parses taxon names in abundance table. Specific for 'new' MG-RAST output 
  ### tables with structure "metagenomeID.DB.taxonlevel.csv"
  ### Table should be read in as tab-separated .csv file.
  ### By Lee Stanish June 27, 2016
  library(stringr)
  
  tax <- x[,2]
  alltax <- str_split_fixed(as.character(tax), pattern=";", 
            n=7)
  colnames(alltax) <- c("domain", "phylum", "class", "order", "family",
                        "genus", "species")
  alltax <- data.frame(alltax)
  return(alltax)
}
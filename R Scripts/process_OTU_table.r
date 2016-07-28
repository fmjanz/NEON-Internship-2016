######### Ordination analysis of otu tables using non-metric MDS #####

#setwd("/Users/stanish/Documents/R")   # for Macs

### Generified script to import an OTU table and prepare for downstream analysis
### LFS June 22 2016

library(vegan)
library(ggplot2)
library(grid)
library("stringr")

### Load required function ###
t.otu <- function(x, s, e) {
  ## Written by Lee Stanish
  ## transpose OTU table for ordination analysis##
  ## x=OTU table; s=starting column; e=ending column
  x <- as.data.frame(x) 
  x1 <- t(x[,c(s:e)])      # transpose OTU counts
  colnames(x1) <- x[,1]     #replace OTU ids as column names 
  x1
}

fill.taxonomy.newMG <- function(x) {
  ### Parses taxon names in abundance table. Specific for 'new' MG-RAST output 
  ### tables with structure "metagenomeID.DB.taxonlevel.csv"
  ### Table should be read in as tab-separated .csv file.
  ### By Lee Stanish June 27, 2016
  
  #tax <- x[,1]
  tax <- x
  alltax <- str_split_fixed(as.character(tax), pattern=";",n=7)
  colnames(alltax) <- c("domain", "phylum", "class", "order", "family","genus", "species")
  alltax <- data.frame(alltax)
  return(alltax)
}
### End required function ##


######    Start of code    ######
setwd("C:/Users/fjanz/Documents/GitHub/NEON-Internship-2016")
#setwd("/Users/stanish/Documents/R")  # set to your dir
otu <- read.csv("Metagenome_beta_div_data.csv",header = TRUE)  # set to your file
env <- read.csv("Metagenome_beta_div_metadata.csv")

## transpose OTU table using function t.otu 
otu.t <- t.otu(x=otu,s=2,e=length(otu))  # adjust numbers as needed for your OTU table
ncol(otu.t)
#rownames(otu.t) <- otu[1,c(2:77)]
taxonomy <- fill.taxonomy.newMG(otu$X)   # check that #columns in otu.t=length of taxonomy


###########           Analysis code (for later)       ####################

###Rarefy OTU table###
### A) Calculate total # seqs for each sample
tototu <- apply(otu.t, MARGIN=1, FUN=sum)
min(tototu)			# rarefy to this value

### B) Rarefy
oturare <- rrarefy(otu.t, min(tototu))

### C) Generate rarefaction curves (test)
rarecurve(otu.t, step=100000,sample=min(tototu))

### D) Remove OTUs with 0 or 1 abundance
otumax <- apply(oturare, 2, max)
otu.1 <- oturare[, otumax >= 2]
taxonomy.1 <- droplevels(taxonomy[otumax >= 2,])

#Transform rarefied OTU table
otu.pct <- decostand(otu.1, method="total")
otu.hell <- decostand(otu.1, method="hellinger")  # Hellinger transformation


######## nMDS  #####
mds <- metaMDS(otu.hell, k=2, distance="bray", trymax=20, autotransform=FALSE)

# Shepard plot to compare actual dissimilarities to computed dissimilarities
diss <- vegdist(otu.hell, distance="bray")
stressplot(mds, diss) ## pretty nice fit

mds.sc <- scores(mds, display="sites", choices=c(1:2))
mds.spec <- scores(mds, display="species", choices=c(1:2))


########     Plot results (not generified)    ########

###All samples#
plot(mds, type="p", choices=c(1,2), cex.lab=1.5,cex.axis=1.5,font=1.5,font.axis=2,font.lab=1.5)

mds.sc <- data.frame(mds.sc)


plot(mds, type="none", choices=c(1,2), cex.lab=1.5,cex.axis=1.5,font=1.5,font.axis=2,
  font.lab=1.5, xlim=c(-0.09,0.12),ylim=c(-0.09,0.09),xlab= "NMDS Axis 1", ylab="NMDS Axis 2",
  main = "Comparison of Sample Types Across Sites", frame.plot=FALSE)
points(mds.sc[,1], mds.sc[,2],display="sites", col=ifelse(env$sampleType=="Individual","darkblue","orangered4"),
       pch = as.integer(env$siteID))
ordihull(mds, env$siteID, lty=2)
legend("topright",pch=c(seq(1,6)), col="black", c("CPER","DSNY","HARV","OSBS","STER","TALL"), 
       bty="o",box.col="black", cex=.8)
legend("bottomright",pch=c(15,15),col=c("darkblue","orangered4"), c("Individual","Composite"),
       bty="o",box.col = "black")

##save to pdf file##
pdf(file="OR_3DnMDS_Taxa10percent.pdf")
plot(mds.sc[,1], mds.sc[,2], pch=20,col="gray", cex=0.8, xlim=c(-1.5,1.5),
     ylim=c(-1,1.3))
text(mdsspec5[,1], mdsspec5[,2], labels=spec5IDs$Class, cex=0.9, 
     col=as.integer(spec5IDs$Phylum))
mtext("Bacterial taxa >= 10% abundance classified to family level")
legend("topleft", bty="n", legend=c(levels(spec5IDs$Phylum)), pt.bg=c(1:7),
       pch=22)
dev.off()



#########  -----------MULTIVARIATE ANALYSES ------------- #########
#Step 1: remove groups with only one sample
to.rm <- c(which(env$siteID=="TALL"))
env.rm <- env[-to.rm,]
otu.hell.rm <- otu.hell[-to.rm,]

diss <- vegdist(otu.hell.rm, method="euclidean")
diss2 <- vegdist(otu.hell.rm, method="bray")
x <- betadisper(diss,group=env.rm$siteID)
x1 <- betadisper(diss2,group=env.rm$siteID)
permutest(x)  
anova(x)
boxplot(x)
TukeyHSD(x)

permutest(x1)  
anova(x1)
boxplot(x1)
TukeyHSD(x1)
#Removing those 2 samples greatly improves stability using either BC or Euc dist

##### PERMANOVA analysis for group membership    #####
adonis(diss~siteID, data=env.rm)
adonis(diss~sampleType, data=env.rm, strata=env.rm$siteID)

adonis(diss2~siteID, data=env.rm)
adonis(diss2~sampleType, data=env.rm, strata=env.rm$siteID)

#########  MANTEL TESTS  #######
envdis <- vegdist(alldat.rm$Free.Cl,method="euclidean")
mantel(envdis,diss,strata=alldat.rm$Municipality)

####  ANOSIM for group membership   ####
mod <- anosim(otu.hell[-52,], grouping=env$ChlType[-52],distance="bray",
              strata=env$State[-52])

#---------------------------------------------------------------------------------------------------------------------

##### Repeat for functional potential #####
otu <- read.csv("Functional_beta_div_data.csv",header=TRUE)
env <- read.csv("Functional_beta_div_metadata.csv", header= TRUE)

otu.t <- t.otu(x=otu,s=2,e=length(otu))  # adjust numbers as needed for your OTU table
ncol(otu.t)

###########           Analysis code (for later)       ####################

###Rarefy OTU table###
### A) Calculate total # seqs for each sample
tototu <- apply(otu.t, MARGIN=1, FUN=sum)
min(tototu)			# rarefy to this value

### B) Rarefy
oturare <- rrarefy(otu.t, min(tototu))

### C) Generate rarefaction curves (test)
rarecurve(otu.t, step=100000,sample=min(tototu))

### D) Remove OTUs with 0 or 1 abundance
otumax <- apply(oturare, 2, max)
otu.1 <- oturare[, otumax >= 2]
#taxonomy.1 <- droplevels(taxonomy[otumax >= 2,])

#Transform rarefied OTU table
otu.pct <- decostand(otu.1, method="total")
otu.hell <- decostand(otu.1, method="hellinger")  # Hellinger transformation


######## nMDS  #####
mds <- metaMDS(otu.pct, k=2, distance="bray", trymax=20, autotransform=FALSE)

# Shepard plot to compare actual dissimilarities to computed dissimilarities
diss <- vegdist(otu.pct, distance="bray")
stressplot(mds, diss) ## pretty nice fit

mds.sc <- scores(mds, display="sites", choices=c(1:2))
mds.spec <- scores(mds, display="species", choices=c(1:2))


########     Plot results (not generified)    ########

###All samples#
plot(mds, type="p", choices=c(1,2), cex.lab=1.5,cex.axis=1.5,font=1.5,font.axis=2,font.lab=1.5)

mds.sc <- data.frame(mds.sc)

#ordination plot
plot(mds, type="none", choices=c(1,2), cex.lab=1.5,cex.axis=1.5,font=1.5,font.axis=2,
     font.lab=1.5, xlim=c(-0.09,0.12),ylim=c(-0.09,0.09),xlab= "NMDS Axis 1", ylab="NMDS Axis 2",
     main = "Collected Functional Gene Composition", frame.plot=FALSE)
points(mds.sc[,1], mds.sc[,2],display="sites", col=ifelse(env$sampleType=="Individual","dodgerblue4","tan3"),
       pch = (as.integer(env$siteID) + 6), bg=ifelse(env$sampleType=="Individual","dodgerblue4","tan3"),
       cex=1.4)
ordihull(mds, env$siteID, lty=2, draw="polygon", col="black", alpha=50)
text(mds,labels= env$siteID)
legend(0.1,0.05,pch=c(seq(7,12)), col="black", c("CPER","DSNY","HARV","OSBS","STER","TALL"), 
       bty="o",box.col="black", cex=1)
legend(0.1,-0.01,pch=c(15,15,4),col=c("dodgerblue4","tan3"), c("Individual","Composite"),
       bty="o",box.col = "black")


######## Select more abundant functions to plot ######
functionMax <- apply(otu.1,2,max)
#hist(functionMax, breaks=24)
mdsfn5000 <- mds$species[functionMax>20000,]   ##nMDS functional coordinates for OTUs
mdsfn5000IDs <- rownames(mdsfn5000)

points(mdsfn5000[,1], mdsfn5000[,2], display="sites", col="dodgerblue4", pch=20,cex=1.6)
text(mdsfn5000[,1], mdsfn5000[,2], display="sites", col="black", labels=seq(1,20), cex=1, pos=3)


##save to pdf file##
pdf(file="OR_3DnMDS_Taxa10percent.pdf")
plot(mds.sc[,1], mds.sc[,2], pch=20,col="gray", cex=0.8, xlim=c(-1.5,1.5),
     ylim=c(-1,1.3))
text(mdsspec5[,1], mdsspec5[,2], labels=spec5IDs$Class, cex=0.9, 
     col=as.integer(spec5IDs$Phylum))
mtext("Bacterial taxa >= 10% abundance classified to family level")
legend("topleft", bty="n", legend=c(levels(spec5IDs$Phylum)), pt.bg=c(1:7),
       pch=22)
dev.off()



#########  -----------MULTIVARIATE ANALYSES ------------- #########
#Step 1: remove groups with only one sample
to.rm <- c(which(env$siteID=="TALL"))
env.rm <- env[-to.rm,]
otu.pct.rm <- otu.pct[-to.rm,]

diss <- vegdist(otu.pct.rm, method="euclidean")
diss2 <- vegdist(otu.pct.rm, method="bray")
x <- betadisper(diss,group=env.rm$siteID)
x1 <- betadisper(diss2,group=env.rm$siteID)
permutest(x)  
anova(x)
boxplot(x)
TukeyHSD(x)

permutest(x1)  
anova(x1)
boxplot(x1)
TukeyHSD(x1)

##### PERMANOVA analysis for group membership    #####
adonis(diss~siteID, data=env.rm)
adonis(diss~sampleType, data=env.rm, strata=env.rm$siteID)

adonis(diss2~siteID, data=env.rm)
adonis(diss2~sampleType, data=env.rm, strata=env.rm$siteID)

#########  MANTEL TESTS  #######
envdis <- vegdist(alldat.rm$Free.Cl,method="euclidean")
mantel(envdis,diss,strata=alldat.rm$Municipality)

####  ANOSIM for group membership   ####
mod <- anosim(otu.hell[-52,], grouping=env$ChlType[-52],distance="bray",
              strata=env$State[-52])




##################################################MLST CODE##############################################


###########################################Set Source for Bioconductor###################################
source("http://bioconductor.org/biocLite.R"); biocLite("BiocUpgrade")
biocLite("Biostrings")
biocLite("seqRFLP")


#################################################1. Load packages #########################################
library(Biostrings)
library(seqRFLP)


#################################### ___LOAD IN DATAFRAMES AND DATABASES___ ############################

###SAMPLE SEQUENCES AND ALLELE SEQUENCES MUST BE ALLIGNED AND TRIMMED TO ENSURE 
###NO DIFFERENCES IN LENGTH, THIS COULD LEAD TO DIFFERENCES BETWEEN SEQEUNCES THAT ARE NOT ACTUALLY 
###SNP's OR INDELS

#### insert check that sequences are same length ####

####################################2. Set WD for each gene ############################################

####rpoB####
rpoB.allele.wd <- "H:\\PhD\\Sequencing\\Fasta\\Sequence Types\\rpoB"
rpoB.seqs.wd <- "H:\\PhD\\Sequencing\\Fasta\\All\\All"
rpoB.results.wd <- "H:\\PhD\\Data\\MLST\\Allele.Numbers\\rpoB"

####ITS####
ITS.allele.wd <- "H:\\PhD\\Sequencing\\Fasta\\Sequence Types\\ITS"
ITS.seqs.wd <- "H:\\PhD\\Sequencing\\Fasta\\All\\All"
ITS.results.wd <- "H:\\PhD\\Data\\MLST\\Allele.Numbers\\ITS"

###################################3. Load allele database for each gene ###############################


####rpoB####
setwd (rpoB.allele.wd) 
rpoB_allele <- readDNAStringSet("Alleles.fas")
rpoB_allele <- rpoB_allele[sort(names(rpoB_allele))]
rpoB_allele.names <- names(rpoB_allele)

rpoB_allele.DF <- data.frame(Allele.No = as.integer(rpoB_allele.names),seq=as.data.frame(rpoB_allele))


####ITS####
setwd (ITS.allele.wd) 
ITS_allele <- readDNAStringSet("Alleles.fas")
ITS_allele <- ITS_allele[sort(names(ITS_allele))]
ITS_allele.names <- names(ITS_allele)

ITS_allele.DF <- data.frame(Allele.No = as.integer(ITS_allele.names),seq=as.data.frame(ITS_allele))

####################################4. Load alignments of sample sequences #############################

####rpoB####

setwd (rpoB.seqs.wd)

rpoB_seqs <- readDNAStringSet("rpoB_all.fas")
rpoB_seqs <- rpoB_seqs[sort(names(rpoB_seqs))]
rpoB_seqs.names <- names(rpoB_seqs)

rpoB_seqs.DF <- data.frame(name=rpoB_seqs.names,seq=as.data.frame(rpoB_seqs), Allele= rep(NA,length(rpoB_seqs.names)))

####ITS####
setwd (ITS.seqs.wd)

ITS_seqs <- readDNAStringSet("ITS_all.fas")
ITS_seqs <- ITS_seqs[sort(names(ITS_seqs))]
ITS_seqs.names <- names(ITS_seqs)

ITS_seqs.DF <- data.frame(name=ITS_seqs.names,seq=as.data.frame(ITS_seqs), Allele= rep(NA,length(ITS_seqs.names)))


###########################################____INSERT____###############################################
###########################################___FUNCTION___###############################################

MLST <- function (gene,allele.wd,alleles,seqs,alleles.DF,seqs.DF,results.wd) 
{

seqs.names <- names(seqs) 

for ( i in 1:length(seqs.names)) {
  
  Allele.match <- alleles[alleles.DF$x==seqs[i],]
  Allele.number <- as.integer(names(Allele.match))
  QC <- length(Allele.number)
  if(QC > 0){
  seqs.DF[i,3] <-Allele.number 
  }
  else{seqs.DF[i,3] <- 0}
  
}
  
##Check for any unknown alleles##
unknown_seqs <- seqs.DF[seqs.DF$Allele==0,] ### find all sequences that weren't assigned
unknown_alleles <- unique(unknown_seqs$x) ### store all unique unknown sequences 
unknowns <- length(unknown_alleles) ### number of unique unkown sequences

message("Number of unknown sequences before update")

print(unknowns) 


## If unknown alleles = 0 then function skips to writing csv of allele profiles ##


##If unkown alleles > 0 then assign allele numbers to unkown alleles##

if(unknowns>0){
  
alleles.DF.new <- alleles.DF  

for(i in 1:length(unknown_alleles))
{
  
 alleles.DF.new[(length(alleles.DF$Allele.No)+i),2] <- unknown_alleles[i]
alleles.DF.new[(length(alleles.DF.new$Allele.No)),1] <- length(alleles.DF.new$Allele.No) 

}


## Overwrite and update fasta file with allele numbers and sequences ##

setwd (allele.wd) ## return to directory of original fasta files

allele.fas <- dataframe2fas(alleles.DF.new) ### Create fasta file of new allele numbers 

write.fasta(allele.fas, "Alleles.fas") ## overwrite previous fasta file to include new allele numbers

alleles.new <- readDNAStringSet("Alleles.fas")
alleles.new <- alleles.new[sort(names(alleles.new))]
alleles.names.new <- names(alleles.new)



for (i in 1:length(seqs.names)) {
  
  Allele.match <- alleles.new[alleles.DF.new$x==seqs[i],]
  Allele.number <- as.integer(names(Allele.match))
  QC <- length(Allele.number)
  if(QC > 0){
    seqs.DF[i,3] <-Allele.number 
  }
  else{seqs.DF[i,3] <- 0}
  }





## Check that there are no more unkown alleles ##
 
  unknown_seqs <- seqs.DF[seqs.DF$Allele==0,] ### find all sequences that weren't assigned
  unknown_alleles <- unique(unknown_seqs$x) ### store all unique unknown sequences 
  unknowns <- length(unknown_alleles) ### number of unique unkown sequences
  
  message("Number of unknown sequences after update")
  
  print(unknowns) 
}


## Write CSV to update database on allele numbers ##

Profiles <- data.frame(Sample=seqs.DF$name,Allele.No=seqs.DF$Allele)

setwd (results.wd)

write.csv(Profiles,"Allele.Profiles.csv")


View(seqs.DF) 
}


#########################################___RUN FUNCTION___############################################

####rpoB####
MLST("rpoB",rpoB.allele.wd,rpoB_allele,rpoB_seqs,rpoB_allele.DF,rpoB_seqs.DF,rpoB.results.wd)
####ITS####
MLST("ITS",ITS.allele.wd,ITS_allele,ITS_seqs,ITS_allele.DF,ITS_seqs.DF,ITS.results.wd)


##########################################################Load Allele Profiles################################################## 

setwd(rpoB.results.wd)
rpoB.profiles <- read.csv("Allele.Profiles.csv",header=TRUE,row.names=1)
#View(rpoB.profiles)

setwd(ITS.results.wd)
ITS.profiles <- read.csv("Allele.Profiles.csv",header=TRUE,row.names=1)
#View(ITS.profiles)

###############################################################################################################################################
#######################################################Assign sequence types###################################################################
###############################################################################################################################################

setwd ("H:\\PhD\\Data\\MLST\\Allele.Numbers")

#1. Sort names of profiles from each loci#################################################################################

rpoB.profiles <- rpoB.profiles[sort(rpoB.profiles$Sample),]
ITS.profiles <- ITS.profiles[sort(ITS.profiles$Sample),]

rpoB.names <- c(as.character(rpoB.profiles$Sample))
ITS.names <- c(as.character(ITS.profiles$Sample))

#2.Quality Control to check that all names match and are in the correct order##################################################################
check1 <- gltA.names==rpoB.names ## checks that sequences are in right order all should return all TRUE if not there is an error 
check2 <- rpoB.names==ITS.names ## checks that sequences are in right order all should return all TRUE if not there is an error  

length(check1[check1==FALSE])
length(check2[check2==FALSE])

## If all checks are 0 then run code below.

#3. Create dataframes containing sequence types for all samples

Allele.profiles <- data.frame(Sample=gltA.names, glta=gltA.profiles$Allele.No, rpoB=rpoB.profiles$Allele.No,ITS=ITS.profiles$Allele.No, ST=(rep(NA,length(ITS.profiles$Allele.No))))

View(Allele.profiles)


Allele.profiles.short <- data.frame(Sample=gltA.names,Profile=rep(NA,length(names)),ST=rep(NA,length(gltA.names)))


for (i in 1:length(gltA.names))
{
  Allele.profiles.short[i,2] <- as.numeric(paste(as.numeric(gltA.profiles[i,2]),as.numeric(rpoB.profiles[i,2]),as.numeric(ITS.profiles[i,2]),sep=""))
}

View(Allele.profiles.short)

#4. Create dataframe of sequence types
setwd ("H:\\PhD\\Data\\MLST\\Allele.Numbers")

ST.Ref <- read.csv("ST.Reference.csv", header=T) 
ST.Ref <- ST.Ref[,-1]


for ( i in 1:length(gltA.names)) {
  
  ST.match <- Allele.profiles.short[i,2] %in% ST.Ref$Profile 
  
  if(ST.match==TRUE){
    ST <-ST.Ref[ST.Ref$Profile==Allele.profiles.short[i,2],] 
    Allele.profiles.short[i,3] <- ST$ST
  }
  else{ Allele.profiles.short[i,3] <- 0}
  
}

unknown.profiles <- Allele.profiles.short[Allele.profiles.short$ST==0,]
unknown.ST <- unique(unknown.profiles$Profile)

length(unknown.ST)


for(i in 1:length(unknown.ST))
{
  
  ST.Ref[(length(ST.Ref$ST)+1),1] <- as.numeric(unknown.ST[i])
  ST.Ref[(length(ST.Ref$ST)),2] <- length(ST.Ref$ST) 
  
}

View(ST.Ref)

for ( i in 1:length(gltA.names)) {
  
  ST.match <- Allele.profiles.short[i,2] %in% ST.Ref$Profile
  
  if(ST.match==TRUE){
    ST <-ST.Ref[ST.Ref$Profile==Allele.profiles.short[i,2],] 
    Allele.profiles.short[i,3] <- ST$ST
  }
  else{ Allele.profiles.short[i,3] <- 0}
  
}

unknown.profiles <- Allele.profiles.short[Allele.profiles.short$ST==0,]
unknown.ST <- unique(unknown.profiles$Profile)

length(unknown.ST)

##########Assign STs to Allele.Profiles################################

for (i in 1:length(Allele.profiles.short$Sample))
{
  Allele.profiles[i,5] <- Allele.profiles.short[Allele.profiles.short$Sample==Allele.profiles[i,1],3]
}

View(Allele.profiles)
## Write output to files ##

setwd ("H:\\PhD\\Data\\MLST\\Allele.Numbers")

write.csv(ST.Ref, "ST.Reference.csv") 

write.csv(Allele.profiles,"ST.Profiles.csv")



 










 














  



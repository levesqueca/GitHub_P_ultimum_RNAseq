###################################################################################################
#Getting started in R
#Set the working directory for getting files
###################################################################################################
setwd("/isilon/biodiversity/users/shared/Pythium_ultimum_RNAseq/")

shared_path_out <- paste(getwd(), "/CAZy_trimmed", sep="")

ref_transcriptomes_folder <- paste(getwd(), "/CAZy_prediction_test/referenceTranscriptomes/", sep="")
###################################################################################################
#Load required library packages
###################################################################################################
source("http://bioconductor.org/biocLite.R")
biocLite("spliceSites")
install.packages("openxlsx")
library(openxlsx)
library(Biostrings)

###################################################################################################
#Getting the fasta seqs for each species:
###################################################################################################

metadata_file2 <- read.table(paste("/isilon/biodiversity/users/shared/Pythium_ultimum_RNAseq/CAZy_prediction_test/referencesMetadataFinal.csv", sep=""),
                         header=TRUE, 
                         sep = ",",
                         skip = 0, 
                         stringsAsFactors = FALSE 
)

metadata_file2 <- metadata_file2[order(metadata_file2$Species),] 

dbCAN_file2 <- read.table(paste(dataPath, "/dbCAN_CAZy_GH131_SLH_allrefs.csv", sep=""),
                                header=TRUE, 
                                sep = ",",
                                skip = 0, 
                                stringsAsFactors = FALSE 
)


i <- 2
#AA_fasta <- AAStringSet()
rm(AA_fasta)
for(i in 1:length(metadata_file2$ProteinsFasta)){
  temp5 = list(readAAStringSet(paste(ref_transcriptomes_folder, metadata_file2$ProteinsFasta[i], sep=""), 
                                            format="fasta"))
  if (i == 1) {
    AA_fasta <- temp5
  } else {
    AA_fasta <- append(AA_fasta,temp5) }
}


attributes(AA_fasta[[1]])

# How can we call up the sequence names???:
for(i in 1:length(AA_fasta)){
  names(AA_fasta[[i]]) = sub(" .*","",names(AA_fasta[[i]]))
}


# Turn CAZy family names into factors in order to create
# a list at the same time, in order to use indexing, use numbering: (cazy and species names)
#dbCAN_file2$family.hmm <- as.factor(dbCAN_file2$family.hmm)

dbCAN_tempFile$family.hmm <- as.factor(dbCAN_tempFile$family.hmm)

species <- as.factor(metadata_file2$Species)

levels(species)
levels(species)[16] # ["Pythium_ultimum_var_ultimum"]

#dbCAN_file2$family.hmm[1]
dbCAN_tempFile$family.hmm[1]

factor(species)

which(species == "Pythium_ultimum_var_ultimum")

i <- 1
j <- 1
dir.create(paste(dataPath, "/trimmedByCAZyDomainAllCAZy", sep = ""), showWarnings = TRUE, recursive = FALSE)
trimmedFastaPath <- paste(dataPath, "/trimmedByCAZyDomainAllCAZy/", sep = "")
for(i in 1:length(levels(dbCAN_tempFile$family.hmm))){
newdata <- subset(dbCAN_tempFile, family.hmm == levels(dbCAN_tempFile$family.hmm)[i] )
rm(x_trim)
  for(j in 1:nrow(newdata)) {
      temp <- AAStringSet(
        AA_fasta[[which(species == newdata$Species[j])]][names(AA_fasta[[which(species == newdata$Species[j])]])==newdata$query.id[j]], 
        start=newdata$query.start[j], end=newdata$query.end[j], width=NA, use.names=TRUE)
      names(temp) <- paste(newdata$Species[j], names(temp), sep="_" )
      if (j == 1) {
        x_trim <- temp
      } else {
        x_trim <- c(x_trim,temp) }
  }
writeXStringSet(x_trim, file=paste(trimmedFastaPath, newdata$family.hmm[1],"_AAtrimmed.fasta", sep=""), append=FALSE, format="fasta") 
}


ID_fasta_files <- list.files(path = trimmedFastaPath, pattern = "\\.fas$|\\.fasta$", recursive = FALSE, full.names = FALSE)


# makes a linsi command from this file without  reorientation
for(i in 1:length(ID_fasta_files)){
  cmd <- paste("/opt/bio/mafft/bin/mafft --reorder --quiet '",  trimmedFastaPath, ID_fasta_files[i], "' > ", 
             trimmedFastaPath, "mafft_aligned_", ID_fasta_files[i], sep = "")
  system(cmd)
}


# The following is to generate a fasta per cazy family that includes sequences
# for these families found inall our reference species. We did this a long time
# ago before re-running the dbCAN cazy pipeline, and we did it only for the oomycetes
# we had from Marcello's paper - so this will refresh and update that for us.
# This way if I want to see the structure of a set of genes from a family of interest,
# I'll be able to look at that and not be limited to only the family domain portions
# when I want to check up on something in these sets.
i <- 1
j <- 1
dir.create(paste(dataPath, "/subsetbyCAZyDomainAllCAZy", sep = ""), showWarnings = TRUE, recursive = FALSE)
parsedFastaPath <- paste(dataPath, "/subsetbyCAZyDomainAllCAZy/", sep = "")
for(i in 1:length(levels(dbCAN_tempFile$family.hmm))){
  newdata <- subset(dbCAN_tempFile, family.hmm == levels(dbCAN_tempFile$family.hmm)[i] )
  rm(x_trim)
  for(j in 1:nrow(newdata)) {
    temp <- AAStringSet(
      AA_fasta[[which(species == newdata$Species[j])]][names(AA_fasta[[which(species == newdata$Species[j])]])==newdata$query.id[j]], 
      width=NA, use.names=TRUE)
    names(temp) <- paste(newdata$Species[j], names(temp), sep="_" )
    if (j == 1) {
      x_trim <- temp
    } else {
      x_trim <- c(x_trim,temp) }
  }
  writeXStringSet(x_trim, file=paste(parsedFastaPath, newdata$family.hmm[1],"_Parsed.fasta", sep=""), append=FALSE, format="fasta") 
}


ID_fasta_files <- list.files(path = parsedFastaPath, pattern = "\\.fas$|\\.fasta$", recursive = FALSE, full.names = FALSE)


# makes a linsi command from this file without  reorientation
for(i in 1:length(ID_fasta_files)){
  cmd <- paste("/opt/bio/mafft/bin/mafft --reorder --quiet '",  trimmedFastaPath, ID_fasta_files[i], "' > ", 
               trimmedFastaPath, "mafft_aligned_", ID_fasta_files[i], sep = "")
  system(cmd)
}


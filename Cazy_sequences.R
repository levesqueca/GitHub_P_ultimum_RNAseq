

library(openxlsx)
library(ape)
library(Biostrings)

remove.packages("xlsx")

shared_path <- "/isilon/biodiversity/users/shared/Pythium_ultimum_RNAseq/"

i <- 3

# Problem because there are two formats of columns, had to do it in two sets 
genome_spp <- c("Pyap","Pyar","Pyir","Pyiw","Pyuu","Pyus","Pyve","Phra","Phso","Phin","Ha")
genome_spp1 <- c("Pyap","Pyus","Phra","Phso","Phin","Ha")
genome_spp2 <- c("Pyar","Pyir","Pyiw","Pyuu","Pyve")

cols_extract1 <- c(1,2,3,15)
cols_extract2 <- c(1,2,3,16)

Cazy_data1 <- data.frame()
for(i in 1:length(genome_spp1)) {
  temp <- read.xlsx(paste(shared_path, "References/Final\ Table\ CAZyme\ 050813.xlsx", sep=""), sheet=genome_spp1[i], 
                    cols=cols_extract1, startRow=2, colNames=TRUE)
  temp$species <-   genome_spp1[i]
  colnames(temp) <- c("CAZy", "Organism", "SequenceID", "Sequence","species")
  Cazy_data1 <-  rbind(temp, Cazy_data1) 
}

Cazy_data2 <- data.frame()
for(i in 1:length(genome_spp2)) {
  temp <- read.xlsx(paste(shared_path, "References/Final\ Table\ CAZyme\ 050813.xlsx", sep=""), sheet=genome_spp2[i], 
                    cols=cols_extract2, startRow=2, colNames=TRUE)
  temp$species <-   genome_spp2[i]
  colnames(temp) <- c("CAZy", "Organism", "SequenceID", "Sequence","species")
  Cazy_data2 <-  rbind(temp, Cazy_data2) 
}

# to put the two sets together
Cazy_data <-  rbind(Cazy_data1, Cazy_data2) 

#To show data with NA
nrow(Cazy_data[is.na(Cazy_data$Sequence),])

Cazy_data_noNA <- na.omit(Cazy_data)

my_AAstring <- AAStringSet(Cazy_data_noNA$Sequence)

my_names <- paste(Cazy_data_noNA$species,Cazy_data_noNA$CAZy, Cazy_data_noNA$SequenceID,sep="_")

names(my_AAstring) <- my_names

length(unique(my_AAstring))

my_AAstring_unique <- my_AAstring[!duplicated(names(my_AAstring))]

my_names2 <- names(my_AAstring_unique)

my_names2 <- gsub("maker-", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("scaffold_", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("contig_", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("fgenesh_", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("masked-", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("fgenesh-", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("gene-", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("snap-", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("snap_", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("abinit-", "",my_names2, ignore.case = FALSE)
my_names2 <- gsub("-mRNA-1", "",my_names2, ignore.case = FALSE)



names(my_AAstring_unique) <- my_names2

writeXStringSet(my_AAstring, file=paste(shared_path,"/CAZy/new_names.fasta",sep=""), append=FALSE, format="fasta") 


CAZy_groups <- unique(Cazy_data_noNA$CAZy)

i <- 1
for(i in 1:length(CAZy_groups)) {
  dir.create(path= paste(shared_path,"/CAZy/",CAZy_groups[i], sep=""), showWarnings = TRUE, recursive = FALSE)
}

i <- 1
for(i in 1:length(CAZy_groups)) {
  temp <-   my_AAstring_unique[grepl(CAZy_groups[i], names(my_AAstring_unique)), ]
  writeXStringSet(temp, file=paste(shared_path,"/CAZy/",CAZy_groups[i],"/",CAZy_groups[i],".fasta",sep=""), append=FALSE, format="fasta") 
}

###############################################
# I did not use this yet

# makes a linsi command from this file without  reorientation
cmd <- paste("/opt/bio/mafft/bin/linsi --reorder '",  ID_Folder, "/", ID_fasta_files, "' > ", ID_Folder, "/", "query_aligned_fasta.fasta", sep = "")

# runs the command on the linux server
system(cmd)


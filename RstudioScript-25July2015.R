############################################################################################################################
#Installing the required packages for R:
# to be done in the R command line and not in R studio
#source("http://www.Bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
############################################################################################################################
#Getting started in R
#set the working directory > setwd("~/")
#check version installed
############################################################################################################################
#Building the STAR Reference Genome Index:
#test1
#test 2
cmd <- paste("mkdir", "GenomeDir")
system(cmd)
STAR_path <- "/home/AAFC-AAC/girouxem/RNASeq/tools/STAR-STAR_2.4.2a/source/STAR"
Pyuu_ref_path <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref_no_mito.fa"
Pyuu_gff3_path <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref_no_mito.gff3"
Pyuu_gtf_path <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref.gtf"
GenomeDir_path <- "/home/AAFC-AAC/girouxem/RNASeq/GenomeDir/"
P <- "Parent"
Gff3_overhang <- 280 #needs to be max len -1, max length of R1 and R2 is ~280? This must be for the pair?
cmd <- paste(STAR_path,"--runMode",
        "genomeGenerate",
        "--genomeDir",GenomeDir_path,
        "--outFileNamePrefix",GenomeDir_path,
        "--genomeFastaFiles",Pyuu_ref_path,
        "--sjdbGTFtagExonParentTranscript",P,
        "--sjdbGTFfile",Pyuu_gff3_path,
        "--sjdbOverhang",Gff3_overhang
        )
cmd
cat(cmd, "\n")
system(cmd)
############################################################################################################################
#Creating the bowtie index for the FASTA reference/transcriptome file:
gff3 <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref_no_mito.gff3"
fa <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref_no_mito.fa"
bowtie2_build_path <- "/opt/bio/bowtie2/bowtie2-build"
cmd <- paste(bowtie2_build_path, " -f ", fa, 
             " /home/AAFC-AAC/girouxem/RNASeq/References","/",paste("Pyuu_ref", sep="_"), sep="")
cmd
system(cmd)
###########################################################################################################################
#Reading the csv file and getting the fasta reads from the illumina miseq folder:
# reads the  csv file with all the info
OTRI <- read.csv("OTRI.csv", stringsAsFactors=FALSE)
# creates a list of file names from the full download path
fs <- basename(OTRI$FastqFilePath)
# to copy from the isilon raw data folder into my working directory, renaming with fs basenames
for(i in 1:nrow(OTRI)) file.copy(OTRI$FastqFilePath[i], fs[i])
# to unzip all files
for(j in 1:length(fs)) {
  cmd <- paste("gunzip", fs[j])
  cat(cmd, "\n")
  system(cmd) # invoke command
}
###########################################################################################################################
#Parsing the Metadata:
parsed_columns <- data.frame(matrix(unlist(strsplit(as.character(OTRI$BaseCallsName), "_|-")), 
                                    nrow=length(OTRI$BaseCallsName), byrow=T), stringsAsFactors = FALSE)
colnames(parsed_columns) <- c("TimePoint","Condition","Library","Illumina_Lane","Read_Direction","Rest")
Metadata <- cbind(parsed_columns[,c(1:3,5)], OTRI)
Metadata$Platform <- "Illumina"
Metadata$ScientificName <- "Pythium ultimum var ultimum"
Metadata$TimePoint <- sub("T","", Metadata$TimePoint, ignore.case = FALSE)
Metadata$RNA_Replicate <- sub("[^0-9$]", "", Metadata$Condition, ignore.case = FALSE)
Metadata$Condition <- sub("^x.*", "Cholesterol", Metadata$Condition, ignore.case = FALSE)
Metadata$Condition <- sub(".*1$|.*2$|.*3$|.*4$", "Control", Metadata$Condition, ignore.case = FALSE)
write.table(Metadata, file = "My_Metadata.csv", append = FALSE, sep =",", col.names=NA)
###########################################################################################################################
#FastqPairedEndValidator
#Validate fastq R1 and R2 pairs are ordered: Make a Metadata table called MetadataRaw that has the raw reads rows collapsed.
library("reshape2")
MetadataRawPairs <- dcast(data = Metadata, LibraryName + Condition + TimePoint + RNA_Replicate ~ Read_Direction, value.var
                    ="BaseCallsName", FUN=c)
MetadataRawPairs$ShortName <- paste(MetadataRawPairs$Condition, MetadataRawPairs$TimePoint, 
                                    MetadataRawPairs$RNA_Replicate, sep=".")
FastqPairedEndValidator_path <- "/home/AAFC-AAC/girouxem/RNASeq/tools/FastqPairedEndValidator.pl"
for(k in 1:nrow(MetadataRawPairs)) {
  cmd = with(MetadataRawPairs, paste(FastqPairedEndValidator_path, R1, R2))}
cmd
sapply(cmd, function(x) system(x))
#Is there a way to save the results /output to a log file so we can have a record the read-pairs were properly ordered?                          
#The output should look something like this:
#total validated mates: 782384 and 782384
#read-pairs are properly ordered
############################################################################################################################
LibraryName <- unique(Metadata$LibraryName, incomparables=FALSE, fromLast=FALSE)
LibraryName
for(j in 1:length(LibraryName)) {
  cmd <- paste("mkdir", LibraryName[j])
  cat(cmd, "\n")
  system(cmd)
}
############################################################################################################################
#SeqPrep: Removing adapter sequences from fastq reads. Do this prior to any other processing to make them easier to detect.
SeqPrep_path <- "/home/AAFC-AAC/girouxem/RNASeq/tools/SeqPrep/SeqPrep"
system(SeqPrep_path)
FwdAdapter <- "ATCTCGTATGCCGTCTTCTGCTTG" 
RevAdapter <- "TAGAGCATACGGCAGAAGACGAAC"
for(k in 1:nrow(MetadataRawPairs)) {
  cmd <- with(MetadataRawPairs, paste(SeqPrep_path, "-f", R1, "-r", R2, "-1", 
    paste("AdapRem_",R1,".gz", sep=""), "-2", paste("AdapRem_",R2,".gz", sep=""),
    "-A", FwdAdapter, "-B", RevAdapter))
}
cmd
sapply(cmd, function(x) system(x))
#Do we want to capture a log file/record of the output results? The output looks like this when done on just 1 pair (T0-2):
# Pairs Processed:  782384
# Pairs Merged:  0
# Pairs With Adapters:	90284
# Pairs Discarded:	6258
# CPU Time Used (Minutes):	5.885167

# To gzip the raw fastq files to manage space - maybe should delete since the raw files are at the original location 
#on Isilon and here it is redundant?
T_raw_zip <- list.files(path = ".", pattern = "^T.*fastq$")
T_raw_zip
for(k in 1:length(T_raw_zip)) {
 cmd <- paste("gzip", T_raw_zip[k])
   cat(cmd, "\n")
   system(cmd)
 }
############################################################################################################################
#Need to unzip the fastq.gz files after adapter removal:
adap_zip <- list.files(path = ".", pattern = "AdapRem_.*fastq.gz$")
adap_zip
for(k in 1:length(adap_zip)) {
  cmd <- paste("gunzip", adap_zip[k])
  cat(cmd, "\n")
  system(cmd)
}
############################################################################################################################
#FastqPairedEndValidator
#Validate fastq R1 and R2 pairs order: Make a Metadata table called MetadataAdapRem that has the raw reads rows collapsed.
for(k in 1:nrow(Metadata)){
  Metadata$Adapters_Removed <- paste("AdapRem",Metadata$BaseCallsName, sep="_")
}
MetadataAdapRem <- dcast(data = Metadata, LibraryName + Condition + TimePoint + RNA_Replicate ~Read_Direction, value.var
                         = "Adapters_Removed", FUN=c)
MetadataAdapRem$ShortName <- paste(MetadataAdapRem$Condition, MetadataAdapRem$TimePoint, MetadataAdapRem$RNA_Replicate, sep=".")
for(k in 1:nrow(MetadataAdapRem)) {
  cmd = with(MetadataAdapRem, paste(FastqPairedEndValidator_path, R1, R2))}
cmd
sapply(cmd, function(x) system(x))
############################################################################################################################
#Processing with PrinSeq:
cmd <- paste("mkdir", "PrinSeq_logs")
system(cmd)
PrinSeq_path <- "perl /home/AAFC-AAC/girouxem/RNASeq/tools/prinseq-lite-0.20.4/prinseq-lite-0.20.4/prinseq-lite.pl" 
nmax <- 1
trim_left <- 10
trim_tail_left <- 5
trim_tail_right <- 5
trim_qual_window <- 3
trim_qual_type <- "mean"
trim_qual_right <- 25 #See about increasing to 30 - see after running with Andre's modified gff.
trim_qual_rule <- "lt"
lc_method <- "dust"
lc_threshold <- 7
out_good <- "processed"
out_bad <- "null"
min_len <- 50
trim_to_len <- 200 #may need to make shorter - need pairs to be ~230 bp because that is the average length of mapped pairs during STAR for our 
#reads. But not sure how to do this - what about using another program to merge pairs? or, try increasing trim_qual_rule to 30.
log <- "processed_log"
for(k in 1:nrow(MetadataAdapRem)) {
  cmd = with(MetadataAdapRem, paste(PrinSeq_path, " -fastq ", R1, " -fastq2 ", R2, 
              " -ns_max_n ", nmax, 
             " -trim_left ", trim_left, 
             " -trim_tail_left ", trim_tail_left, 
             " -trim_tail_right ", trim_tail_right, 
             " -trim_qual_window ", trim_qual_window, 
             " -trim_qual_type ", trim_qual_type, 
             " -trim_qual_right ",  trim_qual_right,
             " -trim_qual_rule ", trim_qual_rule, 
             " -lc_method ", lc_method, 
             " -lc_threshold ", lc_threshold, 
             " -out_good ", paste("Processed",R1, sep="_"), 
             " -out_bad ", out_bad, 
             " -verbose ",  
             " -min_len ", min_len,
             " -trim_to_len ", trim_to_len,
             " -no_qual_header ", 
             " -log ", "/home/AAFC-AAC/girouxem/RNASeq/PrinSeq_logs","/",paste("Processed_log",LibraryName, sep="_"),
             sep=""))
}
cmd
sapply(cmd, function(x) system(x))
list.files(path = "/home/AAFC-AAC/girouxem/RNASeq/PrinSeq_logs/", pattern=glob2rx("Processed_log*"), full.names=T)
ProcessedLogs_PrinSeq <- list.files(path = "/home/AAFC-AAC/girouxem/RNASeq/PrinSeq_logs/", pattern=glob2rx("Processed_log*"), full.names=T)
ProcessedLogs_PrinSeq
for(k in 1:nrow(Metadata)) {
  Metadata$PrinSeq_Processed_fastq_log <- paste("Processed_log",Metadata$LibraryName, sep="_")
}
# glob2rx("Processed*singletons.fastq") - I don't think we want to keep the singletons, correct?
list.files(path=".", pattern=glob2rx("Processed*singletons.fastq"), full.names=T)
processed_singletons <- list.files(path=".", pattern=glob2rx("Processed*singletons.fastq"), full.names=T)
file.remove(processed_singletons)
# to gzip all files - make a list of AdapRem - manage space
list.files(path=".", pattern=glob2rx("AdapRem*.fastq"), full.names=T)
adap_rem_zip <- list.files(path=".", pattern=glob2rx("AdapRem*.fastq"), full.names=T)
adap_rem_zip
for(k in 1:length(adap_rem_zip)) {
   cmd <- paste("gzip", adap_rem_zip[k])
   cat(cmd, "\n")
   system(cmd)
 }
############################################################################################################################
#FastqPairedEndValidator - Running this helps me for now.
#Must rename the Processed fastq due to the name done by PrinSeq:
list.files(path=".", pattern=glob2rx("Processed_AdapRem*_1.fastq"), full.names=T)
ProcessedR1 <- list.files(path=".", pattern=glob2rx("Processed_AdapRem*_1.fastq"), full.names=T)
ProcessedR1
sapply(ProcessedR1, FUN=function(eachPath){
  file.rename(from=eachPath, to=sub(pattern=(".fastq_1.fastq"), replacement=".fastq", eachPath))
})
list.files(path=".", pattern=glob2rx("Processed_AdapRem*_2.fastq"), full.names=T)
ProcessedR2 <- list.files(path=".", pattern=glob2rx("Processed_AdapRem*_2.fastq"), full.names=T)
ProcessedR2
sapply(ProcessedR2, FUN=function(eachPath){
  file.rename(from=eachPath, to=sub(pattern=("_L001_R1_001.fastq_2.fastq"), replacement="_L001_R2_001.fastq", eachPath))
})
list.files(path = ".", pattern = "Processed_.*fastq$")
for(k in 1:nrow(Metadata)){
  Metadata$Processed_Fastq <- paste("Processed_AdapRem",Metadata$BaseCallsName, sep="_")
}
MetadataProcessed <- dcast(data = Metadata, LibraryName + Condition + TimePoint + RNA_Replicate ~ Read_Direction, value.var
                           ="Processed_Fastq", FUN=c)
for(k in 1:nrow(MetadataProcessed)) {
  cmd = with(MetadataProcessed, paste(FastqPairedEndValidator_path, R1, R2))}
cmd
sapply(cmd, function(x) system(x))
############################################################################################################################
#PrinSeq graph fastq reports: Note - this is time-intensive, consider if better to do fastqc with only parameters Andre is
#interested in - or parameters we want for debugging.
#Make a directory where all the graphs will be:
cmd <- paste("mkdir", "PrinSeq_graph_data") 
cat(cmd, "\n")
system(cmd)
#Generate QC graphs for the processed fastq, output the graph files to the graph directory:
for(k in 1:nrow(MetadataProcessed)){
  cmd = with(MetadataProcessed, paste(PrinSeq_path,
        " -fastq ", R1, " -fastq2 ", R2,
        " -verbose ",
        " -out_good ", out_bad,
        " -out_bad ", out_bad,
        " -graph_data ", " /home/AAFC-AAC/girouxem/RNASeq/PrinSeq_graph_data","/",paste(LibraryName,"processed.gd", sep="_"),
        sep=""))
}
cmd
sapply(cmd, function(x) system(x))
list.files(path = "/home/AAFC-AAC/girouxem/RNASeq/PrinSeq_graph_data/", pattern=glob2rx("*_processed.gd"), full.names=T)
Processed_graphs_PrinSeq <- list.files(path = "/home/AAFC-AAC/girouxem/RNASeq/PrinSeq_graph_data/", pattern=glob2rx("*_processed.gd"), full.names=T)
Processed_graphs_PrinSeq
for(k in 1:nrow(Metadata)) {
  Metadata$PrinSeq_Processed_fastq_graph <- paste(Metadata$LibraryName,"processed.gd", sep="_")
}
############################################################################################################################
#FastQC - I'm putting this in because it's simpler than generating all the PrinSeq graphs and having to break to go online,
#Also, the outputs are very similar, and you can do all the most relevant stats.
#Can you take a look at the reads for quality vs length? Maybe after processing I should trim further? Or adjust some of 
#the PrinSeq Parameters?
library("ShortRead")
list.files(path = ".", pattern = "Processed_.*fastq$")
# runs qa with the files in the directory and with the pattern
fqQC = qa(dirPath = ".", pattern = "Processed_.*fastq$", type = "fastq")
report(fqQC, type = "html", dest = "fastqQAreport")
###Performing some stats (optional):
df <- fqQC[["readQualityScore"]]
#[['readCounts']] # inspect read yield
#[['baseCalls']] # inspect base composition
#[['frequentSequences']] # inspect most common sequences
#df <- within(df, acc_sum <- cumsum(density))
df <- within(df, {cumsumDensity <- ave(density, lane, FUN = cumsum)})
df <- within(df, {sumDensity <- ave(density, lane, FUN = sum)})
df$percentile <- df$cumsumDensity / df$sumDensity
library(ggplot2)
pdf(file = "quality_scores_distribution.pdf", width = 7, height =10)
p1 <-   ggplot(df, aes(x=quality, y=density)) +
  geom_line() + 
  facet_wrap( ~ lane, ncol=4) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 6)) 
print(p1)
p2 <-   ggplot(df, aes(x=quality, y=percentile)) +
  geom_line() +
  facet_wrap( ~ lane, ncol=4) +
  theme_bw() +  
  theme(strip.text.x = element_text(size = 6)) 
print(p2)
dev.off()
# checks some stats
ag <- aggregate(density ~ lane, data = df, max)
df.max <- merge(ag, df)
see_cut_off <- subset(df, df$quality > 35 & df$quality < 35.1)
#Looking at the fastQC html, with respect to nucleotide frequency base calls. How do we know that the difference in AT (higher)
#versus GC(lower), is due to actual GC content unique to species and not due to a GC that is introduced during DNA amplification
# by PCR during library preparation - as this is a major source in variation whereby GC-rich regions remain annealed during
#amplification.
###########################################################################################################################
tophat2_path <- "/opt/bio/tophat/bin/tophat2"
# system(tophat2_path) ##Andre - here I can see so many options - surely there's a way to improve the alignment issue we have
gff3 <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref_no_mito.gff3"
bowind = "Pyuu_ref"
bowtie2_path <- "/opt/bio/bowtie2/bowtie2"
cmd = with(MetadataProcessed, paste(tophat2_path, " -G ", gff3,
                             " -p 12 -o ", "/home/AAFC-AAC/girouxem/RNASeq", "/",LibraryName,"/",LibraryName,"_TopHat",
                             " /home/AAFC-AAC/girouxem/RNASeq/References","/",bowind, 
                             " ", R1," ", R2,
                             sep=""))
cmd
#system(cmd)
# I can't get TopHat2 to run in Rstudio - see error:
# [2015-07-25 20:23:51] Beginning TopHat run (v2.0.5)
# -----------------------------------------------
#   [2015-07-25 20:23:51] Checking for Bowtie
# Bowtie 2 not found, checking for older version..
# Error: Bowtie not found on this system.

#made a script as described for making one for bam to sam in previous script from 22july2015.
############################################################################################################################
#Samtools for TopHat ouputs:
#Sort by name, convert to SAM for HT-Seq-count:
cmd = with(MetadataProcessed, paste(samtools1_path, 
                                    " sort", 
                                    " -n ",
                                    "/home/AAFC-AAC/girouxem/RNASeq","/",LibraryName,"/",LibraryName,"_TopHat","/","accepted_hits.bam",
                                    " /home/AAFC-AAC/girouxem/RNASeq","/",LibraryName,"/",LibraryName,"_TopHat","/",LibraryName,"_sn",
                                    sep=""))
cmd                                    
sapply(cmd, function(x) system(x))

cmd = with(MetadataProcessed, paste(samtools1_path, 
                                    " view ", 
                                    " -o ", 
                                    " /home/AAFC-AAC/girouxem/RNASeq","/",LibraryName,"/",LibraryName,"_TopHat","/",LibraryName,"_sn.sam",
                                    " /home/AAFC-AAC/girouxem/RNASeq","/",LibraryName,"/",LibraryName,"_TopHat","/",LibraryName,"_sn.bam",
                                    sep=""))
cmd
sapply(cmd, function(x) system(x))

# Sort BAM for IGV
cmd = with(MetadataProcessed, paste(samtools1_path, 
                                    " sort ", 
                                    " /home/AAFC-AAC/girouxem/RNASeq","/",LibraryName,"/",LibraryName,"_TopHat","/","accepted_hits.bam",
                                    " /home/AAFC-AAC/girouxem/RNASeq","/",LibraryName,"/",LibraryName,"_TopHat","/",LibraryName,"_s",
                                    sep=""))
cmd
sapply(cmd, function(x) system(x))
# Index BAM files for IGV
cmd = with(MetadataProcessed, paste(samtools1_path, 
                                    " index ", 
                                    " /home/AAFC-AAC/girouxem/RNASeq","/",LibraryName,"/",LibraryName,"_TopHat","/",LibraryName,"_s.bam",
                                    sep=""))
cmd
sapply(cmd, function(x) system(x))

#HTSeq-count for TopHat2 hits:
htseq_count_path <- "/home/AAFC-AAC/girouxem/RNASeq/tools/HTSeq-0.6.1/HTSeq-0.6.1/build/scripts-2.7/htseq-count"
system(htseq_count_path)
MetadataProcessed$countf_TopHat = paste(MetadataProcessed$LibraryName, "TopHat2_count", sep=".")
gff3 <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref_no_mito.gff3"
stranded <- "no"
MINAQUAL <- 10
cmd = with(MetadataProcessed, paste(htseq_count_path, 
                                    " -s ", stranded,
                                    " -a ", MINAQUAL,
                                    " --idattr=Parent ",
                                    " /home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"_TopHat","/",LibraryName,"_sn.sam ",
                                    gff3,
                                    " > ",
                                    "/home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"_TopHat","/",MetadataProcessed$countf_TopHat,
                                    sep=""))
cmd
sapply(cmd, function(x) system(x))

library("edgeR")
#Identify the count files and read them into R using readDGE

#Problem here: 
counts_list_TopHat = paste("/home/AAFC-AAC/girouxem/RNASeq/", MetadataProcessed$LibraryName,"/", MetadataProcessed$LibraryName,"_TopHat","/",MetadataProcessed$countf_TopHat, sep="")

counts_list_TopHat
counts_TopHat = readDGE(counts_list_TopHat)$counts
counts_TopHat
 
write.table(counts_TopHat,file=file.path(results,"counts_readDGE_TopHat.annotated.txt"),sep="\t",row.names=T,col.names=T,quote=F) 


#### HTSEQ_COUNT RESULTS 
# 4. Filter weakly expressed and noninformative (e.g., non-aligned) features: 
noint = rownames(counts_TopHat) %in% c("__no_feature","__ambiguous","__too_low_aQual", 
                                "__not_aligned","__alignment_not_unique") 
mean(colSums(counts_TopHat[!noint,])/colSums(counts_TopHat))
## MEAN % of reads map to features
cpms = cpm(counts_TopHat)  ## counts per million
# In edgeR, it is recommended to remove features without  
# at least 1 read per million in n of the samples,  
# where n is the size of the smallest group of replicates,  
keep = rowSums(cpms >1) >=3 & !noint 
dim(counts_TopHat) ## the number of features you started with 
counts_TopHat = counts_TopHat[keep,] 
dim(counts_TopHat) ## count counts of features you have left over after initial filter 
colnames(counts_TopHat) = MetadataProcessed$LibraryName 
#Create a DGEList object (edgeR's container for RNA-seq count data): 
dTopH = DGEList(counts=counts_TopHat, group=MetadataProcessed$Condition) 
#Estimate normalization factors using, RNA composition and adjust for read depth: 
dTopH = calcNormFactors(dTopH)
#Inspect the relationships between samples using a multidimensional scaling (MDS) plot, as shown in Figure 4:
results <- "14-Aug-2015-TopHat-edgeR"
dir.create(results)
pdf(file.path(results,"MDS-edgeR_TopHat.pdf")) 
plotMDS(dTopH, labels=MetadataProcessed$LibraryName, 
        col = rainbow(length(levels(factor(MetadataProcessed$Condition))))[factor(MetadataProcessed$Condition)],cex=0.6, main="MDS") 
dev.off() 


# edgeR - using glm 
#Create a design matrix to specify the factors that are expected to affect expression levels: 
designTopH = model.matrix( ~ Condition, MetadataProcessed) ## samples is your sample sheet, treatment is a column in the sample sheet 
designTopH 


### Here it is pH6 - pH2.5 (so positive fold change indicates higher expression in Normal, negative is higher in Affected) 
#Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood 
d2TopH = estimateGLMCommonDisp(dTopH, designTopH) 
d2TopH = estimateGLMTrendedDisp(d2TopH, designTopH) 
d2TopH = estimateGLMTagwiseDisp(d2TopH, designTopH) 


#plot the mean-variance relationship: 
pdf(file.path(results,"mean.variance-edgeR_TopHat.pdf")) 

# Plot the relationship between mean expression and variance of expression 
plotMeanVar(d2TopH, show.tagwise.vars=TRUE, NBline=TRUE,main="MeanVar") 

# Plot the Biological Coefficient of Variation (as opposed to technical coefficient of variation) 
plotBCV(d2TopH,main="BCV") 
dev.off() 


#Given the design matrix and dispersion estimates, fit a GLM to each feature: 
fTopH = glmFit(d2TopH, designTopH) 

#Perform a likelihood ratio test, specifying the difference of interest 
deTopH = glmLRT(fTopH, coef=2) ## Treatment coefficient 

#Use the topTags function to present a tabular summary of the differential expression statistics 
ttTopH = topTags(deTopH, n=nrow(dTopH)) ## all tags, sorted 
head(ttTopH$table) ## Check result 
table(ttTopH$table$FDR< 0.05) ## the number of "Statistically Differentially Expressed Genes" at an FDR of 0.05 


#Inspect the depth-adjusted reads per million for some of the top differentially expressed genes: 
ncTopH = cpm(dTopH, normalized.lib.sizes=TRUE) 
rnTopH = rownames(ttTopH$table) 
head(ncTopH[rnTopH,order(MetadataProcessed$Condition)],5) 


#Plot the M (log-fold change) versus A (log-average expression) 
degTopH = rnTopH[ttTopH$table$FDR < .05] 
pdf(file.path(results,"smear-edgeR_TopHat.pdf")) 
plotSmear(dTopH, de.tags=degTopH,main="Smear") 
dev.off() 

#Save the result table as a CSV file: 
write.table(ttTopH$table,file=file.path(results,"toptags_edgeR_TopHat.annotated.txt"),sep="\t",row.names=T,col.names=T,quote=F) 

############################################################################################################################
#Using STAR instead of TopHat:
STAR_path <- "/home/AAFC-AAC/girouxem/RNASeq/tools/STAR-STAR_2.4.2a/source/STAR"
GenomeDir_path <- "/home/AAFC-AAC/girouxem/RNASeq/GenomeDir/"
#Running STAR mapping: **I still need to work on setting parameters.
winAnchorMultimapNmax <- 1000
outFilterMultimapNmax <- 1000
outFilterMatchNminOverLread <- 0.4
outFilterScoreMinOverLread <- 0.4
outFilterMismatchNmax <- 100
seedSearchStartLmax <- 15
outFilterScoreMin <- 0
cmd = with(MetadataProcessed, paste(STAR_path, 
              " --genomeDir ", GenomeDir_path,
              " --readFilesIn ", R1, " ", R2,
              " --outFileNamePrefix ", "/home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,
              " --winAnchorMultimapNmax ", winAnchorMultimapNmax,
              " --outFilterMultimapNmax ", outFilterMultimapNmax,
              " --outFilterMatchNminOverLread ", outFilterMatchNminOverLread,
              " --outFilterScoreMinOverLread ", outFilterScoreMinOverLread,
              " --outFilterMismatchNmax ", outFilterMismatchNmax,
              " --seedSearchStartLmax ", seedSearchStartLmax,
              " --outFilterScoreMin ", outFilterScoreMin,
              " --runThreadN 12", sep=""))
cmd
sapply(cmd, function(x) system(x))
#Notes:
#High % unmapped: Too short - rRNA are typically multi-mappers (getplenty in our samples), and if the rRNA repeats are not in the 
#assembly, they will not be mapped and will be reported as "alignment too short".
#Recall that we removed the mitochondrial DNA from our references fasta and gff3 - perhaps only remove the repeats and keep single copy?
#They will map, and we'll know to remove these..?
############################################################################################################################
#Samtools for STAR alignment outputs:
samtools1_path <- "/opt/bio/samtools1/bin/samtools1"
# Convert to BAM for IGV
cmd = with(MetadataProcessed, paste(samtools1_path, 
                                    " view ", 
                                    " -b ",
                                    " -o ", "/home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"Aligned.out.bam",
                                    " /home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"Aligned.out.sam",
                                    sep=""))
cmd
sapply(cmd, function(x) system(x))
# Sort BAM for IGV
cmd = with(MetadataProcessed, paste(samtools1_path, 
                                    " sort ", 
                                    "/home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"Aligned.out.bam",
                                    " /home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"_s",
                                    sep=""))
# Index BAM files for IGV
cmd
sapply(cmd, function(x) system(x))
cmd = with(MetadataProcessed, paste(samtools1_path, 
                                    " index ", 
                                    "/home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"_s.bam",
                                    sep=""))
cmd
sapply(cmd, function(x) system(x))
#IGV run command: $ java -jar /opt/bio/IGV/igv.jar
#make genome
#load _s.bam file
#load gff3 file
############################################################################################################################
#HTSeq-count
htseq_count_path <- "/home/AAFC-AAC/girouxem/RNASeq/tools/HTSeq-0.6.1/HTSeq-0.6.1/build/scripts-2.7/htseq-count"
system(htseq_count_path)
MetadataProcessed$countf = paste(MetadataProcessed$LibraryName, "count", sep=".")
gff3 <- "/home/AAFC-AAC/girouxem/RNASeq/References/Pyuu_ref_no_mito.gff3"
stranded <- "no"
MINAQUAL <- 10
cmd = with(MetadataProcessed, paste(htseq_count_path, 
                                    " -s ", stranded,
                                    " -a ", MINAQUAL,
                                    " --idattr=Parent ",
                                    " /home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", LibraryName,"Aligned.out.sam ",
                                    gff3,
                                    " > ",
                                    "/home/AAFC-AAC/girouxem/RNASeq", "/", LibraryName, "/", MetadataProcessed$countf,
                                    sep=""))
cmd
sapply(cmd, function(x) system(x))
############################################################################################################################
library("edgeR")
#Identify the count files and read them into R using readDGE
counts_list = sapply(file.path("/home/AAFC-AAC/girouxem/RNASeq", MetadataProcessed$LibraryName, sep=""), dir,pattern=".count$",full.names=T)
counts_list
counts = readDGE(counts_list)$counts
counts
#The names of each list is the following 
names(counts_list)
#### HTSEQ_COUNT RESULTS 
# 4. Filter weakly expressed and noninformative (e.g., non-aligned) features: 
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual", 
                                "__not_aligned","__alignment_not_unique") 
mean(colSums(counts[!noint,])/colSums(counts))
## MEAN % of reads map to features
cpms = cpm(counts)  ## counts per million
# In edgeR, it is recommended to remove features without  
# at least 1 read per million in n of the samples,  
# where n is the size of the smallest group of replicates,  
keep = rowSums(cpms >1) >=3 & !noint 
dim(counts) ## the number of features you started with 
counts = counts[keep,] 
dim(counts) ## count counts of features you have left over after initial filter 
colnames(counts) = MetadataProcessed$LibraryName 
#Create a DGEList object (edgeR's container for RNA-seq count data): 
d = DGEList(counts=counts, group=MetadataProcessed$Condition) 
#Estimate normalization factors using, RNA composition and adjust for read depth: 
d = calcNormFactors(d)
#Inspect the relationships between samples using a multidimensional scaling (MDS) plot, as shown in Figure 4:
results <- "14-Aug-2015-Star-edgeR"
dir.create(results)
pdf(file.path(results,"MDS-edgeR.pdf")) 
plotMDS(d, labels=MetadataProcessed$LibraryName, 
         col = rainbow(length(levels(factor(MetadataProcessed$Condition))))[factor(MetadataProcessed$Condition)],cex=0.6, main="MDS") 
dev.off() 


# edgeR - using glm 
#Create a design matrix to specify the factors that are expected to affect expression levels: 
design = model.matrix( ~ Condition, MetadataProcessed) ## samples is your sample sheet, treatment is a column in the sample sheet 
design 


### Here it is pH6 - pH2.5 (so positive fold change indicates higher expression in Normal, negative is higher in Affected) 
#Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood 
d2 = estimateGLMCommonDisp(d, design) 
d2 = estimateGLMTrendedDisp(d2, design) 
d2 = estimateGLMTagwiseDisp(d2, design) 


#plot the mean-variance relationship: 
pdf(file.path(results,"mean.variance-edgeR.pdf")) 

# Plot the relationship between mean expression and variance of expression 
plotMeanVar(d2, show.tagwise.vars=TRUE, NBline=TRUE,main="MeanVar") 

# Plot the Biological Coefficient of Variation (as opposed to technical coefficient of variation) 
plotBCV(d2,main="BCV") 
dev.off() 


#Given the design matrix and dispersion estimates, fit a GLM to each feature: 
f = glmFit(d2, design) 

#Perform a likelihood ratio test, specifying the difference of interest 
de = glmLRT(f, coef=2) ## Treatment coefficient 

#Use the topTags function to present a tabular summary of the differential expression statistics 
tt = topTags(de, n=nrow(d)) ## all tags, sorted 
head(tt$table) ## Check result 
table(tt$table$FDR< 0.05) ## the number of "Statistically Differentially Expressed Genes" at an FDR of 0.05 


#Inspect the depth-adjusted reads per million for some of the top differentially expressed genes: 
nc = cpm(d, normalized.lib.sizes=TRUE) 
rn = rownames(tt$table) 
head(nc[rn,order(MetadataProcessed$Condition)],5) 


#Plot the M (log-fold change) versus A (log-average expression) 
deg = rn[tt$table$FDR < .05] 
pdf(file.path(results,"smear-edgeR.pdf")) 
plotSmear(d, de.tags=deg,main="Smear") 
dev.off() 

#Save the result table as a CSV file: 
write.table(tt$table,file=file.path(results,"toptags_edgeR.annotated.txt"),sep="\t",row.names=T,col.names=T,quote=F) 



### if you have annotation 
# anno <- read.table(annof,sep="\t",header=T,comment.char="",quote="",as.is=T) 
# write.table(data.frame(tt$table,anno[match(rownames(tt$table),anno$Ensembl.Gene.ID),]),file=file.path(results,"toptags_tibia_edgeR.annotated.txt"),sep="\t",row.names=T,col.names=T,quote=F) 
### END OF DIFFERENTIAL EXPRESSION ANALYSIS, ADDITIONAL PROCESSING ON THE TABLE CAN BE PERFORMED 

#Repeat for timepoint comparisons - just curious...:
#da = DGEList(counts=counts, group=MetadataProcessed$TimePoint) 
#da = calcNormFactors(da) 
#pdf(file.path(results,"MDS-timepoint-edgeR.pdf"))
#plotMDS(da, labels=MetadataProcessed$LibraryName, 
#col = rainbow(length(levels(factor(MetadataProcessed$TimePoint))))[factor(MetadataProcessed$TimePoint)],cex=0.6, main="MDS") 
#dev.off() 
#designa = model.matrix( ~ TimePoint, MetadataProcessed) ## samples is your sample sheet, treatment is a column in the sample sheet 
#designa 
#d2a = estimateGLMCommonDisp(da, designa) 
#d2a = estimateGLMTrendedDisp(d2a, designa) 
#d2a = estimateGLMTagwiseDisp(d2a, designa) 
#pdf(file.path(results,"mean.variance-timepoint-edgeR.pdf")) 
#plotMeanVar(d2a, show.tagwise.vars=TRUE, NBline=TRUE,main="MeanVar") 
#plotBCV(d2a,main="BCV") 
#dev.off() 
#fa = glmFit(d2a, design) 
#dea = glmLRT(fa, coef=2) ## Treatment coefficient 
#tta = topTags(dea, n=nrow(d)) ## all tags, sorted 
#head(tta$table) ## Check result 
#table(tta$table$FDR< 0.05) ## the number of "Statistically Differentially Expressed Genes" at an FDR of 0.05 
#nca = cpm(da, normalized.lib.sizes=TRUE) 
#rna = rownames(tta$table) 
#head(nca[rna,order(MetadataProcessed$TimePoint)],5) 
#dega = rna[tta$table$FDR < .05] 
#pdf(file.path(results,"smear-edgeR-timepoint.pdf")) 
#plotSmear(da, de.tags=deg,main="Smear") 
#dev.off() 
#write.table(tta$table,file=file.path(results,"toptags_edgeR.annotated.timepoint.txt"),sep="\t",row.names=T,col.names=T,quote=F) 


# these two commands will show you the first list
counts_list[[1]]
counts_list$samples

# The second list is your counts with $counts name and this why they had the $counts at the end
# here is another way to do this
counts <- counts_list$count
head(counts,5)


#ii Filter weakly expressed and noninformative (e.g., non-aligned) features using a command like:

noint = rownames(counts) %in%
  c("__no_feature", "__ambiguous", "too_low_aQual",
    "__not_aligned", "__alignment_not_unique")

# This is a check to see these values
counts[noint, ]

cpms = cpm(counts)
keep = rowSums(cpms > 0.1) & !noint
#in the part above that has '=0.1' - this is for the number of replicates - I only have 0.1 at the moment for each, but this will change later.
# problem with too many ambiguous reads
counts = counts[keep, ]

#iii Visualize and inspect the count table as follows:
colnames(counts) = Metadata2$shortname
head( counts[,order(Metadata2$Condition)], 5 )

#iv Create a DGEList object (edgeR's container for RNA-seq count data), as follows:
d = DGEList(counts=counts, group=Metadata2$Condition)

# Names of the list in d
names(d)

d$samples
head(d$counts,5)

#v) Estimate normalization factors using:
d = calcNormFactors(d)
# this creates two matrices into a list
d$samples
head(d$counts,5)


#vi) Inspect the relationships between samples using a multidimensional scaling plot, as shown in Figure 4A:
plotMDS(d, labels=Metadata2$shortname,
        col=c("darkgreen","blue")[factor(Metadata2$Condition)])

#vii) Estimate tagwise dispersion (simple design) using:
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

names(d)



#viii) Create a visual representation of the mean-variance relationship using the plotMeanVar (shown in Figure 5A) and plotBCV (Figure 5B) functions, as follows:
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
plotBCV(d)

#ix) Test for differential expression (\classic" edgeR), as follows:
de = exactTest(d, pair=c("Control","Cholesterol"))

# data in list
names(de)
head(de$table,5)
head(de$comparison,5)
head(de$genes,5)



#x) Follow Step 14 B vi)-ix).
#vi) Use the topTags function to present a tabular summary of the differential expression statistics 
#(Note: topTags operates on the output of exactTest or glmLRT, while only the latter is shown here):

tt = topTags(de, n=nrow(d))
names(tt)
head(tt$table)
tt$adjust.method
tt$comparison
tt$test

#vii) Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:

nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
length(rn)
head(nc,5)
head(nc[rn,],5)
head(nc[rn,order(Metadata2$Condition)],5)


#viii) Create a graphical summary, such as an M (log-fold-change) versus A (log-averageexpression) plot, here showing the 
#genes selected as differentially expressed (with a 5% false discovery rate; see Figure 6A):
# we do not have enough data that passed
deg = rn[tt$table$FDR < .9]
length(deg)
plotSmear(d, de.tags=deg)


#ix) Save the result table as a CSV (comma-separated values) le (alternative formats are possible) as follows:
write.csv(tt$table, file="toptags_edgeREG22June2015.csv")

# Function to generate PrinSeq graph files (.gd).
# Options:
# a = Metadata table being referred to for read input
# b = Column name in x with the names of the input fastq files
# c = Prefix - name of directory where graph .gd files are to go
# d = tag to add to the graph file name, gets added as LibraryName.tag.gd
# e = Col name in x with the names of the second read pairs of the input fastq 
#     files if paired end reads are not merged

# Examples:
# i.e., for when you have single reads, R1 and R2:
# cmd <- MakePrinSeqGraphFiles(metadataRawPairs, metadataRawPairs$R1, prefix, "rawGraphs", 
#                              metadataRawPairs$R2)

# i.e, for when you have merged reads:
# cmd <- MakePrinSeqGraphFiles(metadataAdapRM, metadataAdapRM$MergedReads, prefix, "adapRemMerged")

MakePrinSeqGraphFiles <- function(a, b, c, d, e){
  if(missing(e)){
    cmd <- with(a,
                paste(prinSeqPath,
                      " -fastq ", paste(pathFastq, b, sep = ""),
                      " -out_good null ",
                      " -out_bad null ",
                      " -verbose ",
                      " -graph_data ", paste(sharedPathAn, c, "/", LibraryName, 
                                             ".", d, ".gd", sep = ""),
                      sep = ""))
  } else {
    cmd <- with(a,
                paste(prinSeqPath,
                      " -fastq ", paste(pathFastq, b, sep = ""),
                      " -fastq2 ", paste(pathFastq, e, sep = ""),
                      " -out_good null ",
                      " -out_bad null ",
                      " -verbose ",
                      " -graph_data ", paste(sharedPathAn, c, "/", LibraryName, 
                                             ".", d, ".gd", sep = ""),
                      sep = ""))
  }
  return(cmd)
}

# Function is for compressed input when single reads are input only:
MakePrinSeqGraphFiles2 <- function(a, b, c, d){
  cmd = with(a,
             paste("zcat ", pathFastq, b, ".gz", " | ",
                   prinSeqPath,
                   " -fastq stdin ", 
                   " -out_good null ",
                   " -out_bad null ",
                   " -verbose ",
                   " -graph_data ", paste(sharedPathAn, c, "/", LibraryName, ".", d, ".gd", sep=""),
                   sep = ""))
  return(cmd)
}
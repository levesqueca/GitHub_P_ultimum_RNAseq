# # Function for processing single-end, merged and paired-end reads using PrinSeq
#
# # Processing fastq reads with PrinSeq:
# # For single end or merged reads, use PrinSeqProcessSE. Input must be .gz compressed.
# # For paired end reads, use PrinSeqProcessPE. Input must NOT be .gz compressed.

# # Usage:
# # PrinSeqProcessSE(metaRef=<metadata table name>, fastq1=<colname with fastq Name in metaRef>,
# #                  outPutSuffix=<suffix for processed reads>, <options, comma-separated>)
# # PrinSeqProcessPE(metaRef=<metadata table name>, fastq1=<colname with fastq Name in metaRef>,
# #                 fastq2=<colname with R2 fastq reads Name in metaRef>,
# #                 outPutSuffix=<suffix for processed reads>, <options, comma-separated>)

# # Arguments:
# # metaRef = metadata table to refer to for the fastq read names
# # fastq1  = Column name in metaRef with the name of the fastq files If reads are merged, only use
# #           this. Used in PrinSeqProcessSE and PrinSeqProcessPE.
# # fastq2  = Column name in metaRef with the name of the Read 2 fastq files. Used with 
# #           PrinSeqProcessPE for uncompressed fastq files.
# # outPutSuffix = Suffix to be added to the processed reads that are output

# # Optional Arguments:
# # The user can add any of the options they recorded for processing that have been loaded into the
# # global environment. These must be specified AFTER the required options have been provided.

# # The following function will work only with single-end reads, or with merged reads that are .gz
# # compressed. This is because PrinSeq can not handle more than 1 compressed input from stdin.
PrinSeqProcessSE <- function(metaRef, fastq1, outPutSuffix, ...){
  arguments <- list(...)
  argumentList <- paste(arguments, collapse = " ")
  cmd <- with(metaRef,
              paste("zcat ", pathFastq, fastq1, ".gz", " | ", 
                    prinSeqPath, 
                    " -fastq stdin ",
                    argumentList, 
                    " -out_bad ", outBad,
                    " -out_good stdout ",
                    " -no_qual_header ",
                    " -verbose ",
                    " -log ", paste(sharedPathAn, prefix, "/", 
                                    LibraryName, ".", outPutSuffix, ".log", sep = ""),
                    " | gzip > ", paste(pathFastq, LibraryName, ".", 
                                        outPutSuffix, ".fastq.gz", sep = ""),
                    sep = ""))
  return(cmd)
}

# # The following is for paired end reads where the input fastq files are NOT compressed - PrinSeq 
# # cannot take two inputs from stdin and so piping two fastq.gz files is not possible.
PrinSeqProcessPE <- function(metaRef, fastq1, fastq2, outPutSuffix, ...){
  arguments <- list(...)
  argumentList <- paste(arguments, collapse = " ")
  cmd <- with(metaRef,
              paste(prinSeqPath, 
                    " -fastq ",    pathFastq, fastq1,
                    " -fastq2 ",   pathFastq, fastq2,
                    " ", argumentList, 
                    " -out_bad ",  outBad,
                    " -out_good ", paste(pathFastq, LibraryName, ".", outPutSuffix, sep = ""),
                    " -no_qual_header ",
                    " -verbose ",
                    " -log ",      paste(sharedPathAn, prefix, "/", 
                                         LibraryName, ".", outPutSuffix, ".log", sep = ""),
                    sep = ""))
  return(cmd)
}
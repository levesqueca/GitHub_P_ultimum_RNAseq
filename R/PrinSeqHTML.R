# Function to create html files from the PrinSeq graph files - after the name of the .gd files have been recorded in a metadata table:
# a = metadata table
# b = prefix, name of directory where .html files are to go
# c = name of column in a with the .gd filename

MakePrinSeqHTML <- function(a, b, c) {
  cmd <- with(a,
              paste(prinSeqGraphPath,
                    " -verbose ", 
                    " -i ", sharedPathAn, b, "/", c,
                    " -html_all ",
                    sep = ""))
  return(cmd)
}
# i.e.,
#  cmd <- MakePrinSeqHTML(metadataRawPairs, prefix, metadataRawPairs$RawGraphFiles)
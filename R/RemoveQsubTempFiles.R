#Function to remove output files after running qsubs:

RemoveQsubTempFiles <- function(path, prefixSub) {
  system(paste("find ", path, prefixSub, "/", prefixSub, "*", ".sub.", "*", " -delete ", sep = ""))
}
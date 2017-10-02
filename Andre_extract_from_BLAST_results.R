projPath    <- "/isilon/biodiversity/users/shared/Pythium_ultimum_RNAseq/CAZy_prediction_test/"

library(data.table)
ncbiGH15List <- fread(paste(projPath, "GH15_canDB_list_NCBI-nr.txt", sep = ""),
                      sep = "auto", header="auto")
parsedID <- setDT(ncbiGH15List)[, paste0("ID", 1:5) := tstrsplit(ID, "|", type.convert = TRUE, fixed = TRUE)]

giSet <- data.frame(parsedID$ID4)
colnames(giSet) <- "Accession"
length(giSet$Accession)
giSet$Accession <- sub("\\.[1-9]$", "", giSet$Accession)


access.nb <- giSet$Accession



read.GBxml <-
  function (access.nb) 
  {
    

    
    
    library(XML)
    library(RCurl)
    library(ape)
    library(dplyr)
    N <- length(access.nb)
    nrequest <- N%/%400 + as.logical(N%%400)
    # X <- character(0)
    dfX <- data.frame()
    for (i in 1:nrequest) {
      a <- (i - 1) * 400 + 1
      b <- 400 * i
      if (i == nrequest) 
        b <- N
      
      # André changed HTTPs here
      URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                   paste(access.nb[a:b], collapse = ","),
                   "&rettype=gb&retmode=xml", sep = "")
      dfX <- bind_rows(dfX,  xmlToDataFrame(getURL(URL), colClasses = NULL, homogeneous = NA,
                                            collectNames = TRUE, nodes = list(),
                                            stringsAsFactors = default.stringsAsFactors()) ) 
    }
    
    colnames(dfX) <- gsub("\\-","_", colnames(dfX))
    
    
    
    GB_source_feature <- c("source1\\.\\.", "organism", "mol_type", "isolate", "strain", "culture_collection", "specimen_voucher",
                           "clone",  "variety", "mating_type", "host", "country", "db_xreftaxon\\:", "tissue_type",
                           "collection_date", "PCR_primers", "misc_RNA", "collected_by", "isolation_source", 
                           "db_xrefBOLD\\:", "rRNA", "identified_by", "dev_stage", "organelle", "db_xref",
                           "authority","misc_feature1\\.\\.", "lat_lon")
    
    GBSeq_feature_table <-  data.frame(dfX$GBSeq_feature_table, check.names = FALSE,
                                       stringsAsFactors = FALSE)
    
    colnames(GBSeq_feature_table) <- "sources"
    
    
    for (i in 1:length(GB_source_feature)) {
      GBSeq_feature_table$sources <- sub(GB_source_feature[i], paste("\\$",GB_source_feature[i],sep=""),
                                         GBSeq_feature_table$sources)
    }
    
    GBSeq_feature_table$sources <- sub("\\$\\$", "\\$", GBSeq_feature_table$sources)
    
    for (i in 1:length(GB_source_feature)) {
      temp <- sub(paste("^.*\\$", GB_source_feature[i], sep=("")),"", GBSeq_feature_table$sources)
      GBSeq_feature_table[i+1] <- sub("\\$.*$", "", temp)
      colnames(GBSeq_feature_table)[i+1] <- GB_source_feature[i]
    }
    
    GBSeq_feature_table$GBSeq_primary_accession <- sub("^[0-9]*", "", GBSeq_feature_table[,2])
    GBSeq_feature_table$GBSeq_primary_accession <- sub("\\.[0-9]*$", "", GBSeq_feature_table$GBSeq_primary_accession)
    
    
    dfX_features <- merge(dfX[,c("GBSeq_primary_accession","GBSeq_organism","GBSeq_sequence")],
                          GBSeq_feature_table[,c("GBSeq_primary_accession","specimen_voucher","strain","isolate","culture_collection","clone","host","country")], 
                          by = "GBSeq_primary_accession",  all.x = TRUE,  sort = FALSE)
    
    
    dfX_features[dfX_features==""] <- NA
    
    
    for (i in 1:nrow(dfX_features)) {
      if (!is.na(dfX_features$specimen_voucher[i])) {
        dfX_features$isolate_num[i] <- dfX_features$specimen_voucher[i]
      } else if (!is.na(dfX_features$strain[i])) {
        dfX_features$isolate_num[i] <-  dfX_features$strain[i]
      } else if (!is.na(dfX_features$isolate[i])) {
        dfX_features$isolate_num[i] <- dfX_features$isolate[i]
      } else if (!is.na(dfX_features$culture_collection[i])) {
        dfX_features$isolate_num[i] <- dfX_features$culture_collection[i]
      } else if (!is.na(dfX_features$clone[i])) {
        dfX_features$isolate_num[i] <- dfX_features$clone[i]
      } else  {
        dfX_features$isolate_num[i] <- "NA"
      }
    }
    
    
    obj <- vector("list", N)
    for (i in 1:N) {
      obj[[i]] <- unlist(strsplit(dfX_features$GBSeq_sequence[i], NULL))  
    }
    obj <- as.DNAbin(obj)
    names(obj) <- access.nb
    
    
    attr(obj, "accession_num") <- dfX_features$GBSeq_primary_accession
    attr(obj, "species") <- gsub(" ", "_", dfX_features$GBSeq_organism)
    attr(obj, "strain") <- gsub(" ", "_", dfX_features$isolate_num)
    attr(obj, "host") <- gsub(" ", "_", dfX_features$host)
    attr(obj, "country") <- gsub(" ", "_", dfX_features$country)
    
    obj
    #print(names(attributes(obj)))
  }


names(attributes(obj))




RG2 <- read.GB("KJ639177")

write.table(giSet$Accession, file = file.path(paste(projPath, "GH15_AccessionList.txt", sep = "")),
            append = FALSE,
            quote  = FALSE,
            row.names = FALSE,
            col.names = FALSE)



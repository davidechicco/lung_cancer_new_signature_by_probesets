setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("biomaRt", "annotate", "GEOquery")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(listOfBiocPackages)  

# Function that retrieves the annotations of a probeset list.
# Parameters: 
#  probeset_array  vector of affy_hg_u133a probeset ids
#  flank_len       gene flank length (default=25)
#  flank_stream    "downstream" (default) or "upstream"
#  flank_type      "gene" (default), "coding_gene", "transcript", or "coding_transcript"
retrieveAnnotations_GPL11532 <- function(probeset_array, 
                                      flank_len = 25, 
                                      flank_stream = "downstream", 
                                      flank_type = "gene") {
  
  stopifnot(flank_stream %in% c("downstream", "upstream"))
  stopifnot(flank_type %in% c("gene", "coding_gene", "transcript", "coding_transcript"))
  flank_stream <- paste0(flank_stream, "_flank")
  flank_type <- paste0(flank_type, "_flank")
  
  cat("Attention: the function retrieveAnnotations_GPL11532() works only for GPL11532 datasets\n")
  
  # Gene list
  currSpecieMart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = 'hsapiens_gene_ensembl', 
                  host = "mar2016.archive.ensembl.org")
  
  thisAnnotLookup <- 
    getBM(mart = currSpecieMart, 
          attributes = c("affy_hugene_1_0_st_v1", "ensembl_gene_id", "gene_biotype", "external_gene_name", flank_type), 
          filter = c("affy_hugene_1_0_st_v1", flank_stream), 
          values = list(probeset_array, flank_len), 
          uniqueRows = TRUE, 
          checkFilters = FALSE, bmHeader = TRUE)
          
          
  
  return(thisAnnotLookup)
}



theGEOplatform <- "GPL11532"
platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
print("sort(names(thisGEOplatformAnnotations))")
print(sort(names(thisGEOplatformAnnotations)))   


allProbesetsIDs <- thisGEOplatformAnnotations$"ID"

 
annotations <- retrieveAnnotations_GPL11532(allProbesetsIDs)

cat("annotations:\n")
print(annotations[,c(3,1,2)])

probesetID_sequence_symbol <- annotations[,c(3,1,2)]

if(nrow(annotations) < length(probesetsList)) {
    cat("Attention: only ", nrow(annotations), " probes out of ", length(probesetsList), " have annotations\n", sep="")
    probesWithoutAnnotations <-  length(probesetsList) - nrow(annotations)
    cat("and ",probesWithoutAnnotations, " probes do not have annotations\n", sep="")
}



outputFile <- paste0("../data/GPL11532_all_probesets_sequences_symbols_rand", exe_num,".csv")
write.csv(probesetID_sequence_symbol, file=outputFile, row.names=FALSE)
cat("saved file ", outputFile, "\n", sep="") 
setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)

list.of.packages <- c("easypackages", "plyr") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "GEOquery", "biomaRt")

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
retrieveAnnotations_GPL96 <- function(probeset_array, 
                                      flank_len = 25, 
                                      flank_stream = "downstream", 
                                      flank_type = "gene") {
  
  stopifnot(flank_stream %in% c("downstream", "upstream"))
  stopifnot(flank_type %in% c("gene", "coding_gene", "transcript", "coding_transcript"))
  flank_stream <- paste0(flank_stream, "_flank")
  flank_type <- paste0(flank_type, "_flank")
  
  cat("Attention: the function retrieveAnnotations_GPL96() works only for GPL96 datasets\n")
  
  # Gene list
  currSpecieMart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = 'hsapiens_gene_ensembl', 
                  host = "mar2016.archive.ensembl.org")
  
  thisAnnotLookup <- 
    getBM(mart = currSpecieMart, 
          attributes = c("affy_hg_u133a", "ensembl_gene_id", "gene_biotype", "external_gene_name", flank_type), 
          filter = c("affy_hg_u133a", flank_stream), 
          values = list(probeset_array, flank_len), 
          uniqueRows = TRUE, 
          checkFilters = FALSE, bmHeader = TRUE)
          
          
  
  return(thisAnnotLookup)
}

thisGEOplatform <- "GPL96"
platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)

allProbesetsIDs <- (platform_ann_df$"ID")

# cat("ATTENTION: running just in the head() mode\n")
 
annotations <- retrieveAnnotations_GPL96(allProbesetsIDs)

cat("annotations:\n")
print(annotations[,c(3,2,1)])

probesetID_sequence_symbol <- annotations[,c(3,2,1)]

outputFile <- paste0("../results/GPL96_all_probesets_sequences_symbols_rand", exe_num,".csv")
write.csv(probesetID_sequence_symbol, file=outputFile, row.names=FALSE)

if(nrow(annotations) < length(allProbesetsIDs)) {
    cat("Attention: only ", nrow(annotations), " probes out of ", length(allProbesetsIDs), " have annotations\n", sep="")
    probesWithoutAnnotations <-  length(allProbesetsIDs) - nrow(annotations)
    cat("and ",probesWithoutAnnotations, " probes do not have annotations\n", sep="")
}

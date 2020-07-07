setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("biomaRt")

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


# BatyPartialSignature <- c("206592_s_at", "208710_s_at", "210974_s_at", "207988_s_at", "208679_s_at", "213513_x_at", "201975_at", "210716_s_at", "205764_at", "206562_s_at", "208865_at", "208866_at", "208867_s_at", "213086_s_at", "213860_x_at", "205324_s_at", "213937_s_at", "203430_at", "217993_s_at", "218463_s_at", "202959_at", "202960_s_at", "203072_at", "202818_s_at", "202819_s_at", " 213604_at", "203683_s_at", "222230_s_at", "203945_at", "203946_s_at", " 217720_at", "209665_at", "201620_at", "217543_s_at", "219070_s_at", "218202_x_at", "203478_at", "202073_at", "202074_s_at", " 218979_at", "202636_at", "203090_at", "218327_s_at")

# BatyPartialSignature <- c("210974_s_at", "207988_s_at", "208679_s_at", "213513_x_at", "201975_at", "210716_s_at", "205764_at", "206562_s_at", "208865_at", "208866_at", "208867_s_at", "213086_s_at", "213860_x_at", "205324_s_at", "213937_s_at", "203430_at", "217993_s_at", "218463_s_at", "202959_at", "202960_s_at", "203072_at", "202818_s_at", "202819_s_at", " 213604_at", "203683_s_at", "222230_s_at", "203945_at", "203946_s_at", " 217720_at", "209665_at", "201620_at", "217543_s_at", "219070_s_at", "218202_x_at", "203478_at", "202073_at", "202074_s_at", " 218979_at", "202636_at", "203090_at", "218327_s_at")

 
annotations <- retrieveAnnotations_GPL96(c("205764_at"))
cat("annotations:\n")
print(annotations[,1:3])

if(nrow(annotations) < length(BatyPartialSignature)) {
    cat("Attention: only ", nrow(annotations), " probes out of ", length(BatyPartialSignature), " have annotations\n", sep="")
    probesWithoutAnnotations <-  length(BatyPartialSignature) - nrow(annotations)
    cat("and ",probesWithoutAnnotations, " probes do not have annotations\n", sep="")
}
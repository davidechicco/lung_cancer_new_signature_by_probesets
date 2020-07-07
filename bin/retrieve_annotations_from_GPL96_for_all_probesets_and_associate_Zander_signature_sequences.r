setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

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

source("./utils.r")

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

# cat("annotations:\n")
# print(annotations[,c(3,1,2)])

if(nrow(annotations) < length(allProbesetsIDs)) {
    cat("Attention: only ", nrow(annotations), " probes out of ", length(allProbesetsIDs), " have annotations\n", sep="")
    probesWithoutAnnotations <-  length(allProbesetsIDs) - nrow(annotations)
    cat("and ",probesWithoutAnnotations, " probes do not have annotations\n", sep="")
}

ZanderSignatureData <- read.csv("../data/Zander_signature/Zander_signature_GPL6102_probesets_sequences_symbols_rand8176.csv", header = TRUE, sep =",", stringsAsFactors=FALSE)

listOfZanderProbesets <- ZanderSignatureData$"Array_Address_Id"
listOfZanderGeneSymbols <- ZanderSignatureData$"Symbol"
listOfZanderSequences <- ZanderSignatureData$"SEQUENCE"

counter <- 1

nucleotide_overlap_length <- 21

cat("nucleotide_overlap_length: ", nucleotide_overlap_length, "\n", sep="")

totalIterations <- length(listOfZanderSequences) * nrow(annotations)

for(k in 1:length(listOfZanderSequences)){

      thisLongZanderSequence <- listOfZanderSequences[k]

      for(i in 1:nrow(annotations)){
      
	  percCompletion <- counter*100/(nrow(annotations)*length(listOfZanderSequences))
	  if((counter %% round(totalIterations/20))==0) {
	      cat(dec_two(percCompletion), "%  ", sep="")
	  }
	  
	 # cat(thisGPL96Sequence, "\t versus \t", thisLongZanderSequence, "\n",  sep="")

	  this_probeset_ann <- annotations[i,]
	  thisGPL96Sequence <- this_probeset_ann[, c("Flank (Gene)")] 
	  
	  sequenceFound <- grepl(thisGPL96Sequence, thisLongZanderSequence, fixed = TRUE)
	  if(sequenceFound == TRUE) {
	  
	      cat("\n", thisGPL96Sequence, "\t versus \t", thisLongZanderSequence, "\n",  sep="")
	      cat("Found a sequence for the ", listOfZanderGeneSymbols[k], " ", listOfZanderProbesets[k], " gene:\n", sep="")
	      print(this_probeset_ann[, c(3,1,2)])
	      break;
	      
	  } else {
	  
	    thisLongZanderSequenceStartingPart <- substr(thisLongZanderSequence, 1, nucleotide_overlap_length)
	    thisLongZanderSequenceEndingPart <- substr(thisLongZanderSequence, (nchar(thisLongZanderSequence))-nucleotide_overlap_length, nchar(thisLongZanderSequence))
	    
	    thisGPL96SequenceStartingPart <- substr(thisGPL96Sequence, 1, nucleotide_overlap_length)
	    thisGPL96SequenceEndingPart <- substr(thisGPL96Sequence, (nchar(thisGPL96Sequence))-nucleotide_overlap_length, nchar(thisGPL96Sequence))
	    
	    partialSequenceFound <- ( grepl(thisLongZanderSequenceStartingPart, thisGPL96SequenceStartingPart, fixed = TRUE) | grepl(thisLongZanderSequenceEndingPart, thisGPL96SequenceEndingPart, fixed = TRUE) || 
	    grepl(thisLongZanderSequenceStartingPart, thisGPL96SequenceEndingPart, fixed = TRUE) || grepl(thisLongZanderSequenceEndingPart, thisGPL96SequenceStartingPart, fixed = TRUE) )
	    
	    if (partialSequenceFound) {
	    
		  cat("\n", thisGPL96Sequence, "\t versus \t", thisLongZanderSequence, "\n",  sep="")
		  cat("Found a sequence for the ", listOfZanderGeneSymbols[k], " ", listOfZanderProbesets[k], " gene:\n", sep="")
		  print(this_probeset_ann[, c(3,1,2)])
		  break;
	      }
	    
	  }
	  
	  counter <- counter + 1
      }
}
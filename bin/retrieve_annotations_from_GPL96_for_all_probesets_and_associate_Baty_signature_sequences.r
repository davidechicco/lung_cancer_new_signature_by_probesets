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

FTSJ1_sequence <- "AGTGTCCTGATTATGTCCAGAATTTTCCCTAAAGGCAGGGATTCTTAACCTGGATAGAAGCCAGGGGTA"
SLC37A2_sequence <- "TAGCTGCTCTGGTGGGGCCTCGGAGAACCTTTATATCTAGCTCAGGGATTGTAAATACACCCATCGGCAC"

listofBatyGenes <- c("H200003249H200004109", "H200004476", "H200004564", "H200005989", "H200006229", "H200006553", "H200006915", "H200007113", "H200008089", "H200008119", "H200009361", "H200009657", "H200011731", "H200012139", "H200013667", "H200014044", "H200014047", "H200014122", "H200014666", "H200014959", "H200017073", "H200017355", "H200017863", "H200017969", "H300003194", "H300003690", "H300005457", "H300006812", "H300007836", "H300008134", "H300010348", "H300010368", "H300012071", "H300012698", "H300013306", "H300013688", "H300013783", "H300015410", "H300018503", "H300018776", "H300021626", "H300021867", "H300022216")


listofBatySequences <- c("AGTGTCCTGATTATGTCCAGAATTTTCCCTAAAGGCAGGGATTCTTAACCTGGATAGAAGCCAGGGGTA", "TTTCCTTCTCTGGCAGGTATCCCCAGCAGTCAATTAACAATAAGCCAGTATAAAACACCTAAATAACCA", "AGCAGGTACCTCCGGGCGGGCAAGCCTCTGATTTGCACACCTGAAAATCGCGGAATTGAGTTTCGATAG", 
" TGCTTAAGCCTGCAGATACGTAATGTGACCACTGTTTTGTGTTGACAATATTGCTTTATACAGTTCTTG",
 "GACTCCGCTCACTCCCCACGGCCAGCGTGGGCACAGGACTGACCCTTCTTCAGAGATAATGACATTTTA", "AGAGTCTGTTCTACATGGGCCTGCCCTCCTGTGATGGGCAGAGGCTCCTGGTACATCGAGAAGATTCCT", "GACTCAGCAGGGTGACTTGCCTCAGAGGCTATATCCCAGTGGGGGAACAAAGAGGAGCCTGGTAAAAAA",
 "TTCTGGTACCTTCCCCATGGAGGACACTGAAAAGGCTGGGTTGGGGACAGGGAGTATCACTCCATAAGT", "CTTCAGTGCGATCAAAGTTCTACGTGCGAGAGCCGCCGAATGCCAAACCTGACTGGCTGAAAGTTGGGT", "CCTGGCTATACAGTGCATCCTTGAACTGTCAGCCCACAGCAGCAATATGCTTATTCTATCCACATCCCT", "TAGACAACACCAAGTTTACGATGGACTGCGTAGCGCCTACGATCCACGTGTACGAGCACCATGAGAACA",
 "TGTCATGTAGAAGATGCCTAGCCCTGGGGACCTTGTGAAATACGCAGGAACCAGGGATACCATCTGGTC", "CGCAGAGGGATATGTGAAGCCCCAGTCTTGCATTGACATTGTGATTCGCCATGTGGCACCCATTCCCAG", "CCACCCAGGCCTTTAGAGTCAGATGTCTTCATTGAAGATAGAGCCGAAATGACTGTGTTTGTACGGTCT",
 "GAAGCTTCAGCCCTGCACATTTGAACTAGTCACTCTCCCAGACTTGCGTGGGTCAGTTCTTTCTGAGTA", "CACTCTGGTTCACTGCCTCTGTCACTGGTGCAGCCTGGTACCTGGCTGTATTATGCCCTGTCCTCACCA", "GTGACGTCCGGAGGAGGCAGGAAAAGTTTGGAACGGGAGGAGCAGCTGTCCCTGAGAAAATCAAGATCA",
 "CTTAACTCCCTTGGACGGCCAGATATTCTTGTCATGTGTGGAGGGGTGATACCACCTCAGGATTATGAA", "TTGCATTGTGATGTGGTTGGCTGGGGGCCGACATTGTTGCCCTGTTTGCCGGTGGCCTTCTTATAAAAA", "ATTACTTCAAATCCAAACCAGTAGAGACCCCACCTGAACAGAATGGCACCCTCACCTCCCAGCCCAACA", "GAAGAAAAGGAATGTTTCAGCTCCTGAATCAAGACTTACTAGGCAGTCTGGTGGCACCACAGCTTTGCC",
 "GGTGTTCTCTCAATAACCCACCTTTGGAAATGATGTTTGATGTCGGGAAAACTCAACCACCTCTGATGA", "GACTGTTGGAAGCGAAGTGGAAGCACTGAACCTCCAGGTGACATCTCTGTTTAAGGAGCTTCAAGAGGC", 
"GGGCAAGTTCTCCAGCTCCTTCGCCATCGACAGCATCCTGCGCAAGCCCT", "CCATCTAGTTTCTGTGGGGTATGATTATCACCAAGCCAAACATCATCATGGAGTGGAGGTGAAGGAAGT", "CAGTAAAGTAACTGATCCTGAAACACTGAAGAGTTGTGAAACTGTCACTGAAGAACCAAGTCTTCAACAG",
 "CCATGGACCCTCCCTGAGCAACTTCCTGAGCAGGAAGCCGAAGCCCCCAGAGCCATCCTGGCAGCATTGT", "GACAACCGCAGCCACAGCGCCCACTTCTTTGAGTTTCTCACCAAGGAGCTAGCCCTGGGCCAGGACCGGT", "GGATGTGATGCTGGAGAATTATAGGAACCTGGTCTCCCTGGGTGAGGATAACCTCCTGGGGATGTGCCCT",
 "GTTCATCAGCCAAGTGTATGCTCACTGGGGACCACGATGGGCATTCAGCCTGGTGTGTGGAATAATAGTG", "TAGCTGCTCTGGTGGGGCCTCGGAGAACCTTTATATCTAGCTCAGGGATTGTAAATACACCCATCGGCAC", "AAAAGCCACTGGGGGAGGCCAGTCCCGATCAGTCCTTTCTGGCTCCAGCTCCAGATCCCGTGCGGAGCAA", "AAAATGCTAATTTAAAAGGTCAAGTGAAGCTGCTCCTCACGTTTTGGCGTGCCTGCGCTCTCTGCAGGCA",
 "ATCGGTGCTCGGAACTTGCTCAAATCTATAGCAAAGCAGAGAGAAGCTCAACAGCAGCAACTTCAAGCCC", "CAGCGGATGATAATCTCAAGACACCTTCCGAGCCCGATGATAATCTCAAGACACCTTCCGAGCGTCAGCT", " CCGGTTGGCGGGTGTGGCCTGCGGGGTCCCTCTTGGGGTACATGTCTGTTGCTTACTAAG", "CTAAGCATGAATTGAGGAACAGAAGAAGCAGAGCAGATGATCGGAGCAGCATTTGTTTCTCCCCAAATCT", "CTCACTTCCCTTTTCTCCCACCCTCTTTTCCAAGCTGTTTCGCTTTGCAATATATTACTGGTAATGAGTT",
 "GCTTTGGAGGTCCAGTTACCAACTCGCTTACCCCAATTCGGGTACCTGACAGCGCCCGGGGCCACGAGTG",
 "CAGTGGATAACTACGCCGTCCTCCAGCTGGATATAGTGTGGTGTTACTTCCGCCTGGAACAGCTGGAATG", "GTGAGGGTATAAAGTGTAACCATCAGTTAAACCTCTCCTGTCATTCCTGGCTTCCTTGCTTCAGAATTGA",
 "GCATTAAGAGCTGAAACTTGGCTTTTAGCTGCATGGCATGTTAAAGTACCTCCGATGTGGCTGGAAGCTT", 
"CTGGTTCTCCACCCCTCCCCCATACAAAATCCACAACAAAGCGCAGTGGTCTCTTGTGAA", 
"ATTCAGAATCTCCTGTGGGAAAAGAATGTGAGAAACAAGGACAATCACTGCATGGAGGTCATAAGGCTGA")

counter <- 1

nucleotide_overlap_length <- 23

cat("nucleotide_overlap_length: ", nucleotide_overlap_length, "\n", sep="")

for(k in 1:length(listofBatySequences)){

      thisLongBatySequence <- listofBatySequences[k]

      for(i in 1:nrow(annotations)){
      
	  percCompletion <- counter*100/(nrow(annotations)*length(listofBatySequences))
	  if((counter %% 3000)==0) {
	      cat(dec_two(percCompletion), "%  ", sep="")
	  }
	  
	 # cat(thisGPL96Sequence, "\t versus \t", thisLongBatySequence, "\n",  sep="")

	  this_probeset_ann <- annotations[i,]
	  thisGPL96Sequence <- this_probeset_ann[, c("Flank (Gene)")] 
	  
	  sequenceFound <- grepl(thisGPL96Sequence, thisLongBatySequence, fixed = TRUE)
	  if(sequenceFound == TRUE) {
	  
	      cat(thisGPL96Sequence, "\t versus \t", thisLongBatySequence, "\n",  sep="")
	      cat("Found a sequence contained in the ", listofBatyGenes[k], "\n")
	      print(this_probeset_ann[, c(3,1,2)])
	      break;
	      
	  } else {
	  
	    thisLongBatySequenceStartingPart <- substr(thisLongBatySequence, 1, nucleotide_overlap_length)
	    thisLongBatySequenceEndingPart <- substr(thisLongBatySequence, (nchar(thisLongBatySequence))-nucleotide_overlap_length, nchar(thisLongBatySequence))
	    
	    thisGPL96SequenceStartingPart <- substr(thisGPL96Sequence, 1, nucleotide_overlap_length)
	    thisGPL96SequenceEndingPart <- substr(thisGPL96Sequence, (nchar(thisGPL96Sequence))-nucleotide_overlap_length, nchar(thisGPL96Sequence))
	    
	    partialSequenceFound <- ( grepl(thisLongBatySequenceStartingPart, thisGPL96SequenceStartingPart, fixed = TRUE) | grepl(thisLongBatySequenceEndingPart, thisGPL96SequenceEndingPart, fixed = TRUE) || 
	    grepl(thisLongBatySequenceStartingPart, thisGPL96SequenceEndingPart, fixed = TRUE) || grepl(thisLongBatySequenceEndingPart, thisGPL96SequenceStartingPart, fixed = TRUE) )
	    
	    if (partialSequenceFound) {
	    
		  cat(thisGPL96Sequence, "\t versus \t", thisLongBatySequence, "\n",  sep="")
		  cat("Found a PARTIAL sequence contained in the ", listofBatyGenes[k], "\n")
		  print(this_probeset_ann[, c(3,1,2)])
		  break;
	      }
	    
	  }
	  
	  counter <- counter + 1
      }
}
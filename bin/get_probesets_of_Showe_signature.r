setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

list.of.packages <- c("easypackages", "plyr") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("utils.r")

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "GEOquery")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

GSE_code <- "GSE13255"
thisGEOplatform <- "GPL6102"

# gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)
# 
# if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# gene_expression <- as.data.frame(exprs(gset))

# cat("str(gset@phenoData@data)\n")
# print(str(gset@phenoData@data))
# cat("str(gset@phenoData@data)\n")

LABEL_DETECTED <- TRUE
 
 if(LABEL_DETECTED == TRUE) {
 

    # we retrieve the platform details
    platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
    platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
    print("sort(names(platform_ann_df))")
    print(sort(names(platform_ann_df)))

    geneSymbolFlag <- TRUE
    
    emptyGeneSymbol <- ""
    FIRST_GENE_EXPRESSION_INDEX <- 2
    
    
    ourSignatureGeneSymbols <- c("RSF1", "DYRK2", "YY1", "C19orf12", "THEM2", "TRIO", "MYADM", "BAIAP2", "ROGDI", "DNAJB14", "BRE", "TMEM41A", "C9orf64", "FAM110A", "PCNXL2", "REST", "C19orf62", "C13orf27", "ASCC3", "SLC1A5", "PTPLAD1", "MRE11A", "GTPBP10", "N/A", "SERPINI2", "CREB1", "CCDC53", "USP48", "ZSCAN2")

   ourSignatureGeneAccessionNumbers <- c("NM_016578", "NM_003583", "NM_003403", "NM_001031726", "NM_018473", "NM_007118", "NM_001020820", "NM_017450", "NM_024589", "NM_024920", "NM_199191", "NM_080652", "NM_032307", "NM_031424", "NM_014801", "NM_005612", "NM_014173", "NM_138779", "NM_022091", "NM_005628", "NM_016395", "NM_005590", "NM_033107", "BX118737", "NM_006217", "AK126342", "NM_016053", "NM_032236", "NM_001007072")
    
    if(geneSymbolFlag) {

        cat("\n[start] loop for the association of the gene symbols to the probeset ID's: completed \n", sep="")
        # we start from 2 because 1 is the label
        for(k in 1:length(ourSignatureGeneSymbols)) {
        
            currentCompletionPerc <- k*100 / nrow(ourSignatureGeneSymbols)
            kForPrint <- 1000
            if ((k %% kForPrint)==0) { cat(dec_two(currentCompletionPerc), "% ", sep="") }
            

            theseGenes <- platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k] & platform_ann_df$"Accession"==ourSignatureGeneAccessionNumbers[k], ]
            
            if(nrow(theseGenes)==0) {
            
	      theseGenes <- platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k] & platform_ann_df$"Accession"==paste0(ourSignatureGeneAccessionNumbers[k], ".1"), ]
            
            } 
            
            if(nrow(theseGenes)==0) {
            
	      theseGenes <- platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k] & platform_ann_df$"Accession"==paste0(ourSignatureGeneAccessionNumbers[k], ".2"), ]
            
            }
            
            if(nrow(theseGenes)==0) {
            
	      theseGenes <- platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k] & platform_ann_df$"Accession"==paste0(ourSignatureGeneAccessionNumbers[k], ".3"), ]
            
            }
            
            if(nrow(theseGenes)==0) {
            
	      theseGenes <- platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k] & platform_ann_df$"Accession"==paste0(ourSignatureGeneAccessionNumbers[k], ".4"), ]
            
            }
            
            noPrintFlag <- FALSE
            if(nrow(theseGenes)==0) { 
            
		
		# noPrintFlag <- TRUE
		theseGenes <- platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k], ]
		if(nrow(theseGenes)>=1) {
		      cat("Attention: ", ourSignatureGeneSymbols[k], " with accession ",  ourSignatureGeneAccessionNumbers[k] ," is missing. ", sep="" )
		      cat("We selected all the probesets with ", ourSignatureGeneSymbols[k], " symbol, but with different accession number:\n", sep="")
		}
            }
            
            if(nrow(theseGenes)==0) { 
            
		  theseGenes <- platform_ann_df[platform_ann_df$"Accession"==ourSignatureGeneAccessionNumbers[k], ]
		  if(nrow(theseGenes)>=1) {
		      cat("Attention: ", ourSignatureGeneSymbols[k], " with any accession number is missing\n", sep="" )
		      cat("We selected all the probesets with ", ourSignatureGeneAccessionNumbers[k], " accession number, but with different symbols:\n", sep="")
		}
            }
            
             if(nrow(theseGenes)==0) { 
             
	      cat("Attention: ", ourSignatureGeneSymbols[k], " and its accession number are missing\n", sep="" )
	      noPrintFlag <- TRUE
             
             }
            
#             # platform_ann_df[platform_ann_df$"Accession"==ourSignatureGeneAccessionNumbers[k],]
#		   platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k],]
#             # platform_ann_df[platform_ann_df$"Accession"==paste0(ourSignatureGeneAccessionNumbers[k], ".1"),]
#             # platform_ann_df[platform_ann_df$"Accession"==paste0(ourSignatureGeneAccessionNumbers[k], ".2"),]
#             
#             platform_ann_df[platform_ann_df$"Symbol"==ourSignatureGeneSymbols[k] & platform_ann_df$"Accession"==paste0(ourSignatureGeneAccessionNumbers[k], ".1"), ]
            
            
            # print(theseGenes[,c("ID", "Symbol", "Accession", "SEQUENCE")])
            # print(theseGenes[,c( "Symbol", "ID", "Accession")])
            
            if(noPrintFlag == FALSE) {
		# cat(" ", theseGenes$"Symbol", " ", theseGenes$"ID", " ", theseGenes$"Accession", "\n", sep="")
		
		print(theseGenes[,c( "Symbol", "ID", "Accession", "Definition")])
		cat("\n")
            }
            
        }
        cat("\n [end] loop for the association of the gene symbols to the probeset ID's \n ", sep="")
    
    }

    
   


    

}

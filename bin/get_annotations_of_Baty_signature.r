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

thisGEOplatform <- "GPL6650"

# gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)
# 
# if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# gene_expression <- as.data.frame(exprs(gset))

# cat("str(gset@phenoData@data)\n")
# print(str(gset@phenoData@data))
# cat("str(gset@phenoData@data)\n")

    # we retrieve the platform details
platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
print("sort(names(platform_ann_df))")
print(sort(names(platform_ann_df)))


LABEL_DETECTED <- TRUE
 
 if(LABEL_DETECTED == TRUE) {


    geneSymbolFlag <- TRUE
    
    emptyGeneSymbol <- ""
    FIRST_GENE_EXPRESSION_INDEX <- 2
    

    BatySignatureOperonIDs <- c("H300013783", "H200017863", "H200012139", "H200009657", "H200017969", "H200014047", "H200008119", "H300003690", "H200007113", "H300021867", "H200005989", "H300005457", "H200014959", "H300018776", "H300008134", "H200017073", "H300012698", "H200014044", "H200004476", "H200006915", "H300012071", "H300013306", "H300007836", "H200003249", "H200006553", "H200014122", "H200013667", "H200017355", "H300010368", "H200011731", "H300003194", "H300010348", "H300013688", "H300015410", "H200004109", "H200009361", "H200008089", "H200014666", "H300006812", "H300018503", "H200004564", "H300021626", "H200006229", "H300022216")

    
    if(geneSymbolFlag) {

        cat("\n[start] loop for the association of the gene symbols to the probeset ID's: completed \n", sep="")
        # we start from 2 because 1 is the label
        for(k in 1:length(BatySignatureOperonIDs)) {
        
            currentCompletionPerc <- k*100 / nrow(BatySignatureOperonIDs)
            kForPrint <- 1000
            if ((k %% kForPrint)==0) { cat(dec_two(currentCompletionPerc), "% ", sep="") }
            

            theseGenes <- platform_ann_df[platform_ann_df$"ID"==BatySignatureOperonIDs[k], ]
           
            noPrintFlag <- FALSE
             if(nrow(theseGenes)==0) { 
             
	      cat("Attention: ", BatySignatureOperonIDs[k], " probeset is missing\n", sep="" )
	      
             
             }

            if(noPrintFlag == FALSE) {
		
		# print(theseGenes[,c( "ID", "ORF", "SEQUENCE")])
		# print(theseGenes[,c( "ID", "ORF")])
		
		cat(theseGenes[[1]], "; ", sep="")
		cat(theseGenes[[3]], "; ", sep="")
		cat(theseGenes[[6]], ";\n", sep="")
            }
            
        }
        cat("\n [end] loop for the association of the gene symbols to the probeset ID's \n ", sep="")
    
    }	
}

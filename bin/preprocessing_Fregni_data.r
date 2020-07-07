setwd(".")
options(stringsAsFactors = FALSE)

list.of.packages <- c("easypackages") # other packages
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

GSE_code <- "GSE104636"
folderPath <- "../data/data_Fregni/"
thisGEOplatform <- "GPL6244"

gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gene_expression <- as.data.frame(exprs(gset))

print(str(gset))

# # 
# # # we add the labels
# # labels_df_temp <- as.data.frame(gset@phenoData@data$characteristics_ch1.3)
labels_df_temp <- as.data.frame(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1))
 
labels_df <- as.data.frame(t(labels_df_temp))
colnames(labels_df) <- colnames(gene_expression)
rownames(labels_df) <- "lung_cancer" 
gene_expression_with_labels <- rbind(labels_df, gene_expression)
 
# we retrieve the platform details
platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
gene_expression_with_labels$GeneSymbol <- ""
gene_expression_with_labels$ID <- ""
 
gene_expression_with_labels[1,ncol(gene_expression_with_labels)] <- "lung_cancer"
#  
# # # we start from 2 because 1 is the label
# # for(k in 2:nrow(gene_expression_with_labels)) {
# # 
# #     gene_expression_with_labels[k,]$GeneSymbol <- platform_ann_df[platform_ann_df$ID==rownames(gene_expression_with_labels[k,]),]$"gene_assignment"
# # 
# # }
# 
FIRST_GENE_EXPRESSION_INDEX <- 2
SIZE_SPLIT_STRING <- 2
GENE_SYMBOL_INDEX <- 2

cat("\n[start] loop for the association of the gene symbols to the probeset ID's: completed ", sep="")

# we start from 2 because 1 is the label
for(k in FIRST_GENE_EXPRESSION_INDEX:nrow(gene_expression_with_labels)) {
    
    currentCompletionPerc <- k*100 / nrow(gene_expression_with_labels)
    kForPrint <- 1000
    if ((k %% kForPrint)==0) { cat(dec_two(currentCompletionPerc), "% ", sep="") }

    thisAssignment <- NULL
    thisProbesetID <- NULL
    thisProbesetID <- rownames(gene_expression_with_labels[k,])
    # cat("this probeset ID: ", thisProbesetID, "\t", sep="" )
    thisAssignment <- platform_ann_df[platform_ann_df$ID==thisProbesetID, ]$"gene_assignment"

    split_string <- strsplit(thisAssignment, "//")

    if (length(split_string[[1]]) >= SIZE_SPLIT_STRING) {
        thisGeneSymbol_temp <- NULL
        thisGeneSymbol  <- NULL
        thisGeneSymbol_temp <- split_string[[1]][GENE_SYMBOL_INDEX]
        thisGeneSymbol <- gsub(" ", "", thisGeneSymbol_temp, fixed = TRUE)        
      #  cat("\t this gene symbol: ", thisGeneSymbol, "\n", sep="")

       gene_expression_with_labels[k,]$GeneSymbol <- thisGeneSymbol

    } else {

        gene_expression_with_labels[k,]$GeneSymbol <- thisAssignment

    }
}
cat("\n [end] loop for the association of the gene symbols to the probeset ID's ", sep="")


    ourSignature <- c("ADH1A",     "ANK3",      "ARF6",      "BICD2",      "FAM106A",      "FARP2",      "FKBP1A",      "HLA-DQA1",      "N4BP2L1",      "NPAS2",      "P2RX4",     "PCSK5",      "RBP1",      "RPL10L",      "SLC16A4",      "SLC39A14",      "SOX15",      "SPATA20",      "TMEM45A",      "TRIM3",      "TSFM",      "WARS2")

    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$"GeneSymbol" %in% ourSignature, ])
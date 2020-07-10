setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# set.seed(11)

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

targetName <- "lung_cancer"
GSE_code <- "GSE33532"
thisGEOplatform <- "GPL570"
datasetName <-  "Meister"
# folderPath <- paste0("../data/datasets_", thisGEOplatform, "/")
# createFolderCommand <- paste0("mkdir -p ", folderPath, "; \n", sep="")
# cat(createFolderCommand)
# system(createFolderCommand)


gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gene_expression <- as.data.frame(exprs(gset))

cat("str(gset@phenoData@data)\n")
print(str(gset@phenoData@data))
cat("str(gset@phenoData@data)\n")

LABEL_DETECTED <- TRUE
 
 if(LABEL_DETECTED == TRUE) {
 
    # # # # we add the labels
    library("plyr")
    # label_list <- as.numeric(revalue(gset@phenoData@data$"histology:ch1", c("normal lung"=0, "primary lung tumor"=1)))
    
    label_list <- c()
    i <- 1
    for(thisTitle in gset@phenoData@data$"histology:ch1") {
      
	if(grepl("normal lung", thisTitle)) {
	      label_list[[i]] <- 0
	  } else {
	    label_list[[i]] <- 1
	  }
	    i <- i + 1
      }
    
    cat("label_list:\n")
    print(label_list)
    
    labels_df_temp <- as.data.frame(label_list)
    
    labels_df <- as.data.frame(t(labels_df_temp))
    colnames(labels_df) <- colnames(gene_expression)
    rownames(labels_df) <- targetName
    gene_expression_with_labels <- rbind(labels_df, gene_expression)
    
    gene_expression_with_labels$ID <- rownames(gene_expression_with_labels)

#     # we retrieve the platform details
#     platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
#     platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
#     print("sort(names(platform_ann_df))")
#     print(sort(names(platform_ann_df)))

     ourAggregatedSignature <- c("201481_s_at", "203013_at", "201798_s_at", "218741_at", "214261_s_at", "222168_at", "202286_s_at", "219051_x_at", "218380_at", "215066_at", "213243_at", "217724_at", "201738_at", "27152_at", "202418_at", "201729_s_at", "219233_s_at", "219355_at", "212702_s_at", "202627_s_at", "221610_s_at", "213778_x_at", "209985_s_at", "215504_x_at", "214182_at", "215063_x_at", "212110_at", "217715_x_at", "210705_s_at")
    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% ourAggregatedSignature, ])

    print("dim(patients_data_filtered_ourSignature)\n")
    print(dim(patients_data_filtered_ourSignature))

    patients_data_filtered_ourSignature <- rbind(gene_expression_with_labels[1,], patients_data_filtered_ourSignature)

    # outputFileName <- paste0(folderPath,"patients_data_filtered_ourSignature_", exe_num, ".csv")
    # write.table(patients_data_filtered_ourSignature, file=outputFileName, row.names=TRUE, col.names=TRUE, sep=",")
    # cat("saved file  ", outputFileName, "\n")

    patients_data_filtered_ourSignature_t <- as.data.frame(t(patients_data_filtered_ourSignature))

    library("dplyr")
    patients_data_filtered_ourSignature_t <- patients_data_filtered_ourSignature_t %>% dplyr::select(-targetName,targetName)

    patients_data_filtered_ourSignature_t <-  slice(patients_data_filtered_ourSignature_t, 1:(n()-1))
    patients_data_filtered_ourSignature_t[,targetName] <- as.numeric(patients_data_filtered_ourSignature_t[,targetName]) 

    tableOnlyOnesAndZeros <-  rbind(patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$lung_cancer==0,], patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$lung_cancer==1,])

#     outputFileName <- paste0(folderPath, datasetName, "_", GSE_code, "_", thisGEOplatform, "_patients_data_filtered_ourSignature_transpose_", exe_num, ".csv")
#     write.table(tableOnlyOnesAndZeros, file=outputFileName, row.names=FALSE, col.names=TRUE, sep=",")
#     cat("saved file  ", outputFileName, "\n")
    
    patients_data <- tableOnlyOnesAndZeros  
    
    ROSE_OVERSAMPLE <- FALSE
    source("./svm_poly_kernel.r")
}


setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)


# #  Spira
# 
# fileName <- "../data/datasets_GPL96/Spira_GSE4115_GPL96_patients_data_filtered_ourSignature_transpose_3255.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL96/Heiskanen_GSE68465_GPL96_patients_data_filtered_ourSignature_transpose_4501.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/Kohno_GSE31210_GPL570_patients_data_filtered_ourSignature_transpose_6021.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/Kuner_GSE27489_GPL570_patients_data_filtered_ourSignature_transpose_8101.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/Kuo_GSE103888_GPL570_patients_data_filtered_ourSignature_transpose_2002.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/Mitchell_GSE101929_GPL570_patients_data_filtered_ourSignature_transpose_9908.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/Philipsen_GSE19188_GPL570_patients_data_filtered_ourSignature_transpose_2513.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL571/Rotunno_GSE20189_GPL571_patients_data_filtered_ourSignature_transpose_518.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/Rousseaux_GSE30219_GPL570_patients_data_filtered_ourSignature_transpose_9049.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/Tsay_GSE54495_GPL570_patients_data_filtered_ourSignature_transpose_5659.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL96/Wachi_GSE3268_GPL96_patients_data_filtered_ourSignature_transpose_3594.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL96/Want_GSE19027_GPL96_patients_data_filtered_ourSignature_transpose_2490.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- TRUE
fileName <- "../data/datasets_GPL570/Xu_GSE118370_GPL570_patients_data_filtered_ourSignature_transpose_2313.csv"
TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE
# fileName <- "../data/datasets_GPL570/ZhangL_GSE32175_GPL570_patients_data_filtered_ourSignature_transpose_867.csv"
# TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE

targetName <- "lung_cancer"
thisP <- 0.5
execution_number <- 100


cat("fileName: ", fileName, "\n", sep="")
cat("targetName: ", targetName, "\n", sep="")
cat("TRAIN_SET_OVERSAMPLING_SYNTHETIC: ", TRAIN_SET_OVERSAMPLING_SYNTHETIC, "\n", sep="")

list.of.packages <- c("easypackages", "e1071", "PRROC", "dplyr", "pastecs", "class", "gmodels", "kernlab", "ROSE", "klaR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("easypackages")
libraries(list.of.packages)

threshold <- 0.5



source("./confusion_matrix_rates.r")
source("./utils.r")




# checkAtLeastTwoMinorityElements
checkAtLeastTwoMinorityElements <- function(training_data_set, test_data_set, targetName){
    
     original_training_data_set <- training_data_set
     original_test_data_set <- test_data_set
     
     numEleTrain <- nrow(training_data_set)
     numEleTest <- nrow(test_data_set)
    
      labelTrue <- 1
      labelFalse <- 0
      numZerosTrain <-  nrow(training_data_set[training_data_set[,c(targetName)]==labelFalse, ])
      numOnesTrain <- nrow(training_data_set[training_data_set[,c(targetName)]==labelTrue,  ])
      
      numZerosTest <-  nrow(test_data_set[test_data_set[,c(targetName)]==labelFalse, ])
      numOnesTest <- nrow(test_data_set[test_data_set[,c(targetName)]==labelTrue,  ])
      
      cat("numZerosTrain = ", numZerosTrain, "\n", sep="")
      cat("numOnesTrain = ", numOnesTrain, "\n", sep="")
      cat("numEleTrain = ", numEleTrain, "\n", sep="")
      
     cat("\nnumEleTest = ", numEleTest, "\n", sep="")
      cat("numZerosTest = ", numZerosTest, "\n", sep="")
      cat("numOnesTest = ", numOnesTest, "\n", sep="")
      
      # if there are only zeros in the training set, or there are only ones in the training set
      # if(abs(numZerosTrain - numEleTrain) < 2  || abs(numOnesTrain == numEleTrain) < 2 ){
      if(numZerosTrain < 2 ){
            cat(">>> checkAtLeastTwoMinorityElements() activation\n")
      
            # row names of the first 2 elements having label == 1
            
            rowNameFirstOneTest <-  row.names(test_data_set[test_data_set$lung_cancer==0,])[1]
            rowNameSecondOneTest <-  row.names(test_data_set[test_data_set$lung_cancer==0,])[2]
            
            # dataframe with only the 2 elements having label == 1
            rowNameFirstTwoEleTest <-  test_data_set[(row.names(test_data_set)) %in% c(rowNameFirstOneTest, rowNameSecondOneTest), ]
      
            # new test set removing the 2 elements having label == 1
            newTestSet <-  test_data_set[!(row.names(test_data_set)) %in% c(rowNameFirstOneTest, rowNameSecondOneTest), ]
          
            rowNameFirstZeroTrain <-  row.names(training_data_set[training_data_set$lung_cancer==1,])[1]
            
           # cat("rowNameFirstZeroTrain = ", rowNameFirstZeroTrain, "\n", sep="")
            
            rowNameSecondZeroTrain <-  row.names(training_data_set[training_data_set$lung_cancer==1,])[2]
            
            # cat("rowNameSecondZeroTrain = ", rowNameSecondZeroTrain, "\n", sep="")
            
            rowNameFirstTwoEleTrain <-  training_data_set[(row.names(training_data_set)) %in% c(rowNameFirstZeroTrain, rowNameSecondZeroTrain), ]
            
            # new test set removing the 2 elements having label == 0
            newTrainSet <-  training_data_set[!(row.names(training_data_set)) %in% c(rowNameFirstZeroTrain, rowNameSecondZeroTrain), ]          

            newTrainSet <- rbind(newTrainSet, rowNameFirstTwoEleTest)            
            newTestSet <- rbind(newTestSet, rowNameFirstZeroTrain)            
            
#             cat("~~~~ newTrainSet:\n")
#             imbalance_retriever(newTrainSet$lung_cancer)
#             cat("~~~~ newTestSet:\n")
#             imbalance_retriever(newTestSet$lung_cancer)
            
            
            thisResult <- list("newTrainSet" = newTrainSet, "newTestSet" = newTestSet)

      
        }    else {
        
            thisResult <- list("newTrainSet" = original_training_data_set, "newTestSet" = original_test_data_set)
        }
        
            return(thisResult)
    }


# data read
patients_data <- read.csv(file=fileName,head=TRUE,sep=",",stringsAsFactors=FALSE)


# let's put the target label last on the right 
patients_data <- patients_data%>%dplyr::select(-targetName,targetName)

NUM_METRICS <- 9
confMatDataFrame <- matrix(ncol=NUM_METRICS, nrow=1)
colnames(confMatDataFrame) <-  c("MCC", "F1_score", "accuracy", "TP_rate", "TN_rate", "PPV", "NPV", "PR_AUC", "ROC_AUC")

# info about the dataset
dataset_dim_retriever(patients_data)
imbalance_retriever(patients_data[,c(targetName)])



cat("Number of executions = ", execution_number, "\n", sep="")
# for(exe_i in 1:execution_number)
exe_i <- 1
while(exe_i <= execution_number)
{
    cat("(exe_i = ", exe_i,")\n", sep="")

    patients_data <- patients_data[sample(nrow(patients_data)),] # shuffle the rows

    totalElements <- dim(patients_data)[1]

    subsets_size <- totalElements

    if (subsets_size != totalElements) {
        cat("!!! ATTENTION: We are running the method on a subset of the original dataset, \n", sep="")
        cat("!!! containing only ", subsets_size, " elements \n", sep="")
        cat("!!! instead of ", totalElements, " elements \n", sep="")
    }

    patients_data <- patients_data[1:subsets_size, ]

    dataset_dim_retriever(patients_data)
    imbalance_retriever(patients_data[,c(targetName)])


    target_index <- dim(patients_data)[2]

    training_set_perce <- 70
    cat("training_set_perce = ", training_set_perce, "% \n", sep="")
    test_set_perce <- 100 - training_set_perce
    cat("test_set_perce = ", test_set_perce, "% \n", sep="")

    # the training set is the first 60% of the whole dataset
    training_set_first_index <- 1 # NEW
    training_set_last_index <- round(dim(patients_data)[1]*training_set_perce/100) # NEW
    
    # the test set is the last 20% of the whole dataset
    test_set_first_index <- round(dim(patients_data)[1]*(training_set_perce)/100)+1 # NEW
    test_set_last_index <- dim(patients_data)[1] # NEW    
    
    patient_data_ones <- patients_data[patients_data[,c(targetName)]==1, ]
    patient_data_zeros <- patients_data[patients_data[,c(targetName)]==0, ]
    
    
    cat("[Creating the training set and test set]\n")   
    patients_data_train_complete <- patients_data[training_set_first_index:training_set_last_index, 1:(target_index)] # NEW
    patients_data_test_complete <- patients_data[test_set_first_index:test_set_last_index, 1:(target_index)] # NEW
    
    balancerOutput <- checkAtLeastTwoMinorityElements(patients_data_train_complete,  patients_data_test_complete, targetName)
    patients_data_train_complete <- balancerOutput$newTrainSet
    patients_data_test_complete <- balancerOutput$newTestSet
    
    
    
     if(TRAIN_SET_OVERSAMPLING_SYNTHETIC == TRUE)
         {
            
            cat("[ROSE] We use ROSE oversampling with thisP=", thisP, "\n", sep="")
            # formula
            allFeaturesFormula <- as.formula(paste(as.factor(colnames(patients_data)[target_index]), '.', sep=' ~ ' ))
                  
            data.rose <- ROSE(allFeaturesFormula, data = patients_data_train_complete, p=thisP, seed = 1)$data
            patients_data_train_complete <- data.rose
         }         
    

    cat("[Creating the subsets for the values]\n")
    patients_data_train <- patients_data_train_complete[, 1:(target_index-1)] # NEW
    patients_data_test <- patients_data_test_complete[, 1:(target_index-1)] # NEW


    cat("[Creating the subsets for the labels \"1\"-\"0\"]\n")
    patients_data_train_labels <- patients_data_train_complete[, target_index] # NEW
    patients_data_test_labels <- patients_data_test_complete[, target_index]   # NEW

    # info about the dataset
    cat("### imbalance training set: \n")
    imba_train_flag <- imbalance_retriever(patients_data_train_labels)
    
        # info about the dataset
    cat("### imbalance test set: \n")
    imba_test_flag <- imbalance_retriever(patients_data_test_labels) 

    
    if(imba_train_flag & imba_test_flag) {

            # c_array = c(0.001, 0.01, 0.1, 1, 10)
            fixed_c <- 0.1        
            fixed_exp <- 2

            cat("\n[Training the SVM model (with the hyper-parameter C=",fixed_c,") on training set & applying the SVM to the test set]\n", sep="")
            #patients_data_test_pred <- knn(train = patients_data_train, test = patients_data_test, cl = patients_data_train_labels, k=bestK)

            # svm_model_new <- svm(patients_data_train_labels ~ ., cost=fixed_c, data=patients_data_train, method = "C-classification", kernel = "linear")
            # svm_model_new <- svm(patients_data_train_labels ~ ., cost=fixed_c, data=patients_data_train, method = "C-classification", kernel = "polynomial", degree=fixed_exp)
            svm_model_new <- svmlight(patients_data_train_labels ~ ., data = patients_data_train, svm.options = "-c 0.1 -t 1 -d 2 -r 0")

	    cat("\n\n\nSVM training completed\n\n\n")
            
            patients_data_test_pred <- predict(svm_model_new, patients_data_test)
            
            patients_data_test_binary_predictions <- as.numeric((patients_data_test_pred$posterior[,2])>=threshold)
            
            thisConfMat <- confusion_matrix_rates(patients_data_test_labels, patients_data_test_binary_predictions, "@@@ Test set @@@")
            
            if (exe_i == 1)  confMatDataFrame <-  thisConfMat
            else  confMatDataFrame <- rbind(confMatDataFrame, thisConfMat)

    } else {
    
            cat("# Skip this execution of non-binary condition #\n")
            exe_i <- exe_i - 1
   
    }    
    
    
    exe_i <- exe_i + 1
}

 cat("\n\n\n=== final results ===\n")
 
 cat("Number of executions = ", nrow(confMatDataFrame), "\n", sep="")
 # statistics on the dataframe of confusion matrices
 statDescConfMatr <- stat.desc(confMatDataFrame)
meanRowResults <- (statDescConfMatr)[c("mean"),]
cat("\n\n")
print(dec_three(statDescConfMatr))
cat("\n\n")
print(dec_three(meanRowResults))
cat("\n\n=== === === ===\n")

computeExecutionTime()

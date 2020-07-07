#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
set -o nounset -o pipefail -o errexit
set -o xtrace


i=1
outputFile=""

random_number=$(shuf -i1-10000 -n1)
outputFolder="../results_other_signatures/derived_Sheng_signature"
mkdir -p $outputFolder
outputFile=$outputFolder"/derived_Sheng_signature_rand"$random_number".txt"

Rscript preparation_and_classification_Heiskanen_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Kohno_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Kuner_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Kuo_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Mitchell_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Philipsen_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Rotunno_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Rousseaux_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Spira_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Tsay_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Wachi_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Want_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_Xu_data_probesets.r >> $outputFile; 
Rscript preparation_and_classification_ZhangL_data_probesets.r >> $outputFile; 

grep "^mean" $outputFile

echo -e "\nThe end\n"


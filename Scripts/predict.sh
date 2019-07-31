#!/bin/bash

#purpose: Run REGAd3p+(58 features)
#author: Sumit Tarafder

input_path="../Input/list";
source_code_path="../Codes";


chmod a+x run_RBSURpred.sh

for file in `cat $input_path/id_list.txt` ;
do
    cd ../Scripts
    ./run_RBSURpred.sh $file
    cd $source_code_path
    fname="predictasa$1"
    g++ -o $fname predictASA_processOutput.cpp
   
    ./$fname $file 58 1 1 1
    
    echo "Prediction Complete for ID: $file ."
done
rm $fname

cd ../Scripts



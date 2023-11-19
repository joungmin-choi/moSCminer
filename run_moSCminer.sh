#!/bin/bash

tar -zxvf training_data.tar.gz
tar -zxvf testing_data.tar.gz

# Preprocess Gene expression dataset
python3 preprocessing_gene_expression.py

python3 get_promoter_region.py
./run_bedtools.sh

python3 preprocessing_accessibility.py
python3 preprocessing_step2.py 
python3 merge_omics.py


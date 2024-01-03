# moSCminer: a cell subtype classification framework based on the attention neural network integrating the single-cell multi-omics dataset on the cloud
moSCminer is a web-based platform with cloud computing for cell subtype identification using the single-cell multi-omics datasets. For multi-omics integration, each omics dataset will be converted to a gene-based matrix, considering the biological relationships between each omics. Attention-based neural network model is constructed, where the self-attention module will be applied to each preprocessed omics dataset to transform each feature to a new representation considering the relative importance weight. Representations from each omics will then be concatenated and delivered to the classification module to predict the subtype of each cell. 

The web platform can be accessed from [this link](http://203.252.206.118:5568/)

User can also download the source codes and run the moSCminer model in ther local computer following the below instructions.

## Requirements
* Python (>= 2.7)
* Tensorflow (v1.8.0)
* Bedtools (v2.26.0)
* Python packages : numpy, pandas, os, sys, scikit-learn

### Installation
To install the above requirments, please run below commands in the terminal :
```
pip install tensorflow==1.8.0 # For CPU
pip install tensorflow-gpu==1.8.0 # For GPU
pip install numpy pandas os sys
pip install -U scikit-learn
```

## Usage
Clone the repository or download source code files.

## Inputs
[Note!] All these inputs should be located in the cloned repository.

### 1) Training and Testing dataset (training_data.tar.gz, testing_data.tar.gz)
All omics dataset needs to be compressed into one "tar.gz" formatted file, and each training and testing dataset should be named as **training_data.tar.gz** and **testing_data.tar.gz**.

Please follow the belowing format for each omics dataset :

* Gene expression Dataset
  -  Should be one matrix formatted *.csv file with raw read count values.
  -  Row: Sample, Column: Feature (The first column needs to be 'sampleName')

* DNA methylation, DNA accessibility Dataset
  -  Each sample dataset should be in bed-formatted file. (*.bed,*.wig)

To help your understanding, the example dataset can be downloaded from [here](https://drive.google.com/drive/folders/1OmaDmttI7Ki0M0pgIQBxggszYsNu1zOg?usp=sharing)

### 2) Metadata for training and testing dataset (training_metadata.csv, testing_metadata.csv)
Metadata examples can be found in the repository.

* Training metadata format
  - Column names: OmicsType,FileName,SampleName,CellType
  - For gene expression, leave the "SampleName" and "CellType" column empty.

* Testing metadata format
  - Column names: OmicsType,FileName,SampleName
  - For gene expression, leave the "SampleName" column empty.

## Run moSCminer
Run the below command :
```
./run_moSCminer.sh
```

After running the moSCminer, **results** directory will be created and you will get below outputs :
* /results/test_data_celltype_prediction_results.csv
  - Cell type prediction results for testing dataset
    
* /results/gene_importance_rank.csv
  - Ranking for each gene based on the learned relative importance values from moSCminer
    
* /results/cpg_cluster_importance_rank.csv
  - Ranking for each CpG cluster based on the learned relative importance values from moSCminer
    
* /results/accessibility_importance_rank.csv
  - Ranking for each chromatin accessibility based on the learned relative importance values from moSCminer

## moSCminer (omics-extensible version)
In our paper, the moSCminer was tested using single-cell multi-omics datasets, which comprised two or three distinct omics profiles. However, it is essential to note that our proposed method is not limited to only three omics types. moSCminer allows its application to multi-omics datasets with varying numbers of omics data. This flexibility is achieved by omics-level attention, by employing a self-attention module for each omics dataset, which generates new representations based on the learned relative importance of features. Subsequently, these features are concatenated and used as input for cell subtype classification.

To use the extensible version of moSCminer,
1. Edit the top line 4-6 of **moSCminer_extensible_ver.py** with your paramters
```
vi moSCminer_extensible_ver.py

'''
num_omics = 3  (Number of your omics)
dirname = "./example_dataset/" (directory path having your input files)
resDir = "./results/"  (directory path to save results)
'''
```

2. Prepare the input files following the below rules :
* train dataset filenames should be as following "final_1_train.csv", "final_2_train.csv", "final_3_train.csv" ...
* test dataset filenames should be as following "final_1_test.csv", "final_2_test.csv", "final_3_test.csv" ...
* label files for train and test dataset should be "final_celltype_onehot_train.csv" and "final_celltype_onehot_test.csv", respectively.
* Dataset files 
  -  Should be one matrix formatted *.csv file with raw read count values.
  -  Row: Sample, Column: Feature
* Label files
  -  Should be one-hot encoded matrix
  -  Row: Sample, Column: Cell subtype
 
3. Run the below command
```
python moSCminer_extensible_ver.py
```

All the results will be saved in the *resDir* path you specified.

## Contact
If you have any questions or problems, please contact to **joungmin AT vt.edu**.

import numpy as np
import pandas as pd
import sys
import os

dirname = './'
train_dirname = './training_data/'
test_dirname = "./testing_data/"     
meta_train = pd.read_csv(dirname + 'training_metadata.csv')
meta_test = pd.read_csv(dirname + 'testing_metadata.csv')

train_gene_fileName = meta_train[meta_train['OmicsType'] == 'gene expression']['FileName'].tolist()[0]
test_gene_fileName = meta_test[meta_test['OmicsType'] == 'gene expression']['FileName'].tolist()[0]
train_file = train_dirname + train_gene_fileName
test_file = test_dirname + test_gene_fileName

train_raw = pd.read_csv(train_file, index_col = 0)     
test_raw = pd.read_csv(test_file, index_col = 0)  

# select common samples
meta_train_common_sample = meta_train[meta_train.duplicated(subset='SampleName')]['SampleName']
meta_test_common_sample = meta_test[meta_test.duplicated(subset='SampleName')]['SampleName']
meta_train_common_sample = set(meta_train_common_sample) & set(train_raw.index)
meta_test_common_sample = set(meta_test_common_sample) & set(test_raw.index)

train = train_raw.loc[meta_train_common_sample]
test = test_raw.loc[meta_test_common_sample]

# remove genes which have 0 value for all cells 
def remove_zero_gene(data):
    rm_gene= data[data.sum(axis=0)[data.sum(axis=0)> 0].index]
    return rm_gene

# library size normalization + log transformation
def normalize(data):
    normalize_data = data.apply(lambda x: np.log((x*10000.0/sum(x))+1.0), axis=1)
    return normalize_data

train_rm_gene = remove_zero_gene(train)     
test_rm_gene = remove_zero_gene(test)

train_norm = normalize(train_rm_gene)
test_norm = normalize(test_rm_gene)

# select common features between train and test dataset
train_feature_list = train_norm.columns
test_feature_list = test_norm.columns

common_feature = set(train_feature_list) & set(test_feature_list)

train_data = train_norm[common_feature]
test_data = test_norm[common_feature]

dirsave = './preprocessing_rna/'
try :
    os.mkdir(dirsave)
except :
    print("Exist!")
train_data.to_csv(dirsave+'gene_expression_train.csv', sep=',') 
test_data.to_csv(dirsave+'gene_expression_test.csv', sep=',') 
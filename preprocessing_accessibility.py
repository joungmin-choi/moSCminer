import numpy as np
import pandas as pd
import sys
import os

acc_dir = "./preprocessing_acc/"

# sum gene region for each file
dir_filelist_train = './training_acc_filelist.txt'
dir_filelist_test = './testing_acc_filelist.txt'

dirname_train = acc_dir + 'results_find_common_region/'
dirname_test = acc_dir + 'results_find_common_region/'

dirsave_train = acc_dir + 'results_sum_gene_region_acc_train/'
dirsave_test = acc_dir + 'results_sum_gene_region_acc_test/'

try :
    os.mkdir(dirsave_train)
    os.mkdir(dirsave_test)
except :
    print("Exist!")

dirmeta = './' 
meta_train = pd.read_csv(dirmeta + 'training_metadata.csv')
meta_train = meta_train[meta_train['OmicsType'] == 'DNA accessibility']
meta_test = pd.read_csv(dirmeta + 'testing_metadata.csv')
meta_test = meta_test[meta_test['OmicsType'] == 'DNA accessibility']

filelist_train= pd.read_csv(dir_filelist_train, header=None)[0]
filelist_test= pd.read_csv(dir_filelist_test, header=None)[0]

def sum_gene_region(filelist, dirsave, dirname):
    for file in filelist:
        print(file)
        raw = pd.read_csv(dirname + file + '.bed', sep='\t', header=None)

        groups = raw.groupby([4, 7])
        dna_sum = groups[3].sum()   
        dna_sum_df = pd.DataFrame(dna_sum)
        dna_sum_df.reset_index(inplace=True)
        dna_sum_df.columns = ['chrom', 'gene', 'DNA_accessibility']

        # meta information (region)
        meta = raw[[0, 1, 2, 3, 7]]
        meta.columns = ['chrom', 'Gene start (bp)', 'Gene end (bp)', 'DNA_accessibility', 'gene']

        dna_sum_df.to_csv(dirsave+file+'.csv', sep=',', index=False)
        meta.to_csv(dirsave+file+'_meta.csv', sep=',', index=False)

# Save files        
sum_gene_region(filelist_train, dirsave_train, dirname_train)
sum_gene_region(filelist_test, dirsave_test, dirname_test)

    
# interaction gene region for all file    
dirname_train = dirsave_train
dirname_test = dirsave_test

def interact_gene_region(filelist, dirname):
    filename = filelist[0]
    data = pd.read_csv(dirname + filename+'.csv')
    data = pd.DataFrame(data.set_index(['chrom', 'gene'])['DNA_accessibility'])
    data.columns = [filename]
    
    for filename in filelist[1:]:
        data_next = pd.read_csv(dirname + filename+'.csv')
        data_next = pd.DataFrame(data_next.set_index(['chrom', 'gene'])['DNA_accessibility'])
        data_next.columns = [filename]
        data = pd.concat([data, data_next], join='inner', axis=1)   # intersection
        print('file name : '+filename +' common feature : ' + str(len(data)))
        
    return data
    
data_train = interact_gene_region(filelist_train, dirname_train)
data_test = interact_gene_region(filelist_test, dirname_test)
data_train.columns = meta_train[meta_train['FileName'].isin(data_train.columns)]['SampleName']
data_test.columns = meta_test[meta_test['FileName'].isin(data_test.columns)]['SampleName']

# select common features between train and test dataset
common_features = set(data_train.index) & set(data_test.index)

data_train.reset_index(['chrom', 'gene'], inplace = True)
data_test.reset_index(['chrom', 'gene'], inplace = True)

common_features = pd.DataFrame(common_features)
common_features = common_features.rename(columns = {0: 'chrom', 1: 'gene'})

data_train_common = pd.merge(common_features, data_train, on = ['chrom', 'gene'])
data_test_common = pd.merge(common_features, data_test, on = ['chrom', 'gene'])

#data_train_common = data_train.loc[common_features]
#data_test_common = data_test.loc[common_features]

data_train_common.to_csv(dirsave_train+'intersection_train.csv', mode = "w", index = False)
data_test_common.to_csv(dirsave_test+'intersection_test.csv', mode = "w", index = False)


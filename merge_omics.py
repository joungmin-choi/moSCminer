import pandas as pd
import os 
from sklearn.preprocessing import MinMaxScaler

methyl_dir = "./preprocessing_methyl/"
acc_dir = "./preprocessing_acc/"
gene_dir = "./preprocessing_rna/"

resDir = './final_preprocessing/'
try :
	os.mkdir(resDir)
except :
	print('exist!')

dataTypeList = ['train', 'test']
for dataType in dataTypeList :
	methyl = pd.read_csv(methyl_dir + "methyl_filtered_gene_avg_promoter_" + dataType + ".csv", index_col = 0)
	gene = pd.read_csv(gene_dir + "gene_expression_" + dataType + ".csv", index_col = 0)
	acc = pd.read_csv(acc_dir + "results_sum_gene_region_acc_" + dataType + "/intersection_" + dataType + ".csv")
	del acc['chrom']
	acc.set_index('gene', inplace = True, drop = True)
	acc = acc.T
	if dataType == "train" :
		sample_info = pd.read_csv("./" + dataType + "ing_metadata.csv")
		sample_info = sample_info.dropna(axis = 0)
		celltype_info = sample_info[['SampleName', 'CellType']]
		celltype_info = celltype_info.drop_duplicates(['SampleName'], keep = 'first')
		celltype_info.set_index('SampleName', inplace = True, drop = True)
		gene_feature_df = pd.DataFrame({'feature' : gene.columns.tolist()})
		methyl_feature_df = pd.DataFrame({'feature' : methyl.columns.tolist()})
		acc_feature_df = pd.DataFrame({'feature' : acc.columns.tolist()})
		final_feature_df = pd.merge(methyl_feature_df, acc_feature_df)
		final_feature_df = pd.merge(final_feature_df, gene_feature_df)
		final_feature_df = final_feature_df.drop_duplicates(['feature'], keep = 'first')
	
	gene = gene.T
	gene_filtered = pd.merge(gene, final_feature_df, left_index = True, right_on = "feature")
	gene_filtered = gene_filtered.drop_duplicates(['feature'], keep = 'first')
	gene_filtered.set_index('feature', inplace = True, drop = True)
	gene_filtered = gene_filtered.T 

	methyl = methyl.T 
	methyl_filterd = pd.merge(methyl, final_feature_df, left_index = True, right_on = "feature")
	methyl_filterd= methyl_filterd.drop_duplicates(['feature'], keep = 'first')
	methyl_filterd.set_index('feature', inplace = True, drop = True)
	methyl_filterd.index = methyl_filterd.index + '_cpg'
	methyl_filterd = methyl_filterd.T 

	acc = acc.T 
	acc_filtered = pd.merge(acc, final_feature_df, left_index = True, right_on = "feature")
	acc_filtered = acc_filtered.drop_duplicates(['feature'], keep = 'first')
	acc_filtered.set_index('feature', inplace = True, drop = True)
	acc_filtered.index = acc_filtered.index + "_acc"
	acc_filtered = acc_filtered.T

	scaler = MinMaxScaler()
	gene_filtered_normalized = scaler.fit_transform(gene_filtered)
	gene_filtered_normalized = pd.DataFrame(gene_filtered_normalized)
	gene_filtered_normalized.columns = gene_filtered.columns 
	gene_filtered_normalized.index = gene_filtered.index

	scaler = MinMaxScaler()
	acc_filtered_normalized = scaler.fit_transform(acc_filtered)
	acc_filtered_normalized = pd.DataFrame(acc_filtered_normalized)
	acc_filtered_normalized.columns = acc_filtered.columns
	acc_filtered_normalized.index = acc_filtered.index

	final = pd.merge(gene_filtered_normalized, methyl_filterd, left_index = True, right_index = True)
	final = pd.merge(final, acc_filtered_normalized, left_index = True, right_index = True)

	if dataType == "train" :
		final = pd.merge(final, celltype_info, left_index = True, right_index = True)
		final.to_csv(resDir + "final_multiomics_dataset_with_celltype_" + dataType + ".csv", mode = 'w', index = True)

	final_gene = final.T[:len(gene_filtered_normalized.columns)].T
	final_gene.to_csv(resDir + "final_gene_" + dataType + ".csv", mode = "w", index = True)

	final_methyl = final.T[len(gene_filtered_normalized.columns):len(gene_filtered_normalized.columns)+len(methyl_filterd.columns)].T
	final_methyl.to_csv(resDir + "final_methyl_" + dataType + ".csv", mode = "w", index = True)

	if dataType == "train" :
		final_acc = final.T[len(gene_filtered_normalized.columns)+len(methyl_filterd.columns):-1].T
	else :
		final_acc = final.T[len(gene_filtered_normalized.columns)+len(methyl_filterd.columns):].T
	final_acc.to_csv(resDir + "final_acc_" + dataType + ".csv", mode = "w", index = True)

	if dataType == "train" :
		final_celltype = pd.get_dummies(final['CellType'])
		final_celltype.to_csv(resDir + "final_celltype_onehot_" + dataType + ".csv", mode = "w", index = True)
	else :
		tmp = final_celltype[:1].reset_index(drop = True)
		fake_celltype_df = pd.DataFrame()
		for i in range(len(final_acc)) :
			fake_celltype_df = pd.concat([fake_celltype_df, tmp], axis = 0)
		fake_celltype_df.reset_index(inplace = True, drop = True)
		fake_celltype_df.to_csv(resDir + "final_celltype_onehot_test.csv", mode = "w")



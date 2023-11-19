import os 
import pandas as pd
from sklearn.impute import SimpleImputer
import numpy as np

methyl_dir = "./preprocessing_methyl/"
acc_dir = "./preprocessing_acc/"
gene_dir = "./preprocessing_rna/"

train_gene = pd.read_csv(gene_dir + "gene_expression_train.csv", index_col = 0)
rna_gene_list = train_gene.columns.tolist()
del(train_gene)

train_acc = pd.read_csv(acc_dir + "results_sum_gene_region_acc_train/intersection_train.csv")
acc_gene_list = train_acc['gene'].tolist()
del(train_acc)

common_genes = set(rna_gene_list) & set(acc_gene_list)
rna_gene_list = list(common_genes)

gene_count_dict = {}
for i in range(len(rna_gene_list)) :
	gene_count_dict[rna_gene_list[i]] = 0

filelist = pd.read_csv("training_methyl_filelist.txt", header = None)
dirname = methyl_dir + "results_find_common_region/"
for i in range(len(filelist)) : 
	filename = dirname + filelist[0][i] + ".bed"
	data = pd.read_csv(filename, header = None, sep = '\t')
	bed_gene_list = data[7].unique().tolist()
	for gene in bed_gene_list :
		try :
			gene_count_dict[gene] += 1
		except :
			continue

gene_list = []
count_list = []
for key,val in gene_count_dict.items() :
	gene_list.append(key)
	count_list.append(val)

gene_count_df = pd.DataFrame({'gene' : gene_list, 'count' : count_list})
gene_count_df = gene_count_df.sort_values(by='count', ascending = False)
gene_count_df_filtered = gene_count_df[gene_count_df['count'] > len(filelist)*0.8]
gene_count_df_filtered_gene_list = gene_count_df_filtered['gene'].tolist()

sample_info = pd.read_csv("./training_metadata.csv")
sample_info = sample_info.dropna(axis = 0)

sample_info_unique = sample_info[['SampleName', 'CellType']]
sample_info_unique = sample_info_unique.drop_duplicates(['SampleName'], keep = 'first')
cell_group = sample_info_unique.groupby('CellType').count()
celltype = cell_group.index.tolist()
num_celltype = cell_group['SampleName'].tolist()

celltype_dict = {}
for cell in celltype :
	celltype_dict[cell] = []
	for gene in gene_count_df_filtered_gene_list :
		celltype_dict[cell].append(0)

gene_idx_list = []
for gene_idx in range(len(gene_count_df_filtered_gene_list)) :
	gene_idx_list.append(gene_idx)

gene_count_df_filtered_gene_df = pd.DataFrame({'gene' : gene_count_df_filtered_gene_list, 'gene_idx' : gene_idx_list})

for i in range(len(filelist)) : 
	filename = dirname + filelist[0][i] + ".bed"
	data = pd.read_csv(filename, header = None, sep = '\t')
	data_gene_list = data[7].unique().tolist() 
	tmp_sample_info = sample_info[sample_info['FileName'] == filelist[0][i]]
	sample_name = tmp_sample_info['SampleName'].tolist()[0]
	data_celltype = tmp_sample_info['CellType'].tolist()[0]
	data_gene_df = pd.DataFrame({'gene' : data_gene_list})
	tmp_df = pd.merge(data_gene_df, gene_count_df_filtered_gene_df, on = "gene")
	for j in range(len(tmp_df)) :
		celltype_dict[data_celltype][tmp_df['gene_idx'][j]] += 1

celltype_df = pd.DataFrame(celltype_dict)
celltype_df['gene'] = gene_count_df_filtered_gene_list
#celltype_df.to_csv("celltype_gene_count.csv", mode = "w", index = False)

filtered_cell = []
for i in range(len(celltype_df)) :
	flag_rm = 0
	for j in range(len(celltype)) :
		if celltype_df[celltype[j]][i] < 0.95 * num_celltype[j] : #Maybe need to control the criteria
			flag_rm = 1
			break
	if flag_rm == 0 :
		filtered_cell.append(celltype_df['gene'][i])

filtered_cell_df = pd.DataFrame({'gene' : filtered_cell})
celltype_df_20p = pd.merge(celltype_df, filtered_cell_df, left_on = "gene", right_on = "gene")

resdir = methyl_dir + "filtered_cpg_info/"

try :
	os.mkdir(resdir)
except :
	print("exist!")

final_grouped_gene_df = pd.DataFrame()
for i in range(len(filelist)) : #len(filelist) 
	filename = dirname + filelist[0][i] + ".bed"
	data = pd.read_csv(filename, header = None, sep = '\t')
	del data[4]
	del data[5]
	del data[6]
	data = data.rename(columns = {0 : 'chr', 1 : 'start', 2 :'end', 3 : filelist[0][i], 7 : 'gene'})
	data_filterd = pd.merge(filtered_cell_df, data, how = 'left', on = 'gene')
	data_filterd.to_csv(resdir + filelist[0][i] + ".csv", mode = "w", index = False)
	del data_filterd['chr']
	del data_filterd['start']
	del data_filterd['end']
	grouped_gene = data_filterd.groupby('gene')
	grouped_gene_df = pd.DataFrame(grouped_gene.mean())
	if i == 0 :
		final_grouped_gene_df = grouped_gene_df.copy()
	else :
		final_grouped_gene_df = pd.merge(final_grouped_gene_df, grouped_gene_df, left_index = True, right_index = True)

final_grouped_gene_df = final_grouped_gene_df.T
tmp = final_grouped_gene_df.copy()

final_grouped_gene_df = pd.merge(final_grouped_gene_df, sample_info, left_index = True, right_on = "FileName")

del final_grouped_gene_df['OmicsType']
del final_grouped_gene_df['FileName']
final_grouped_gene_df.set_index("SampleName", inplace = True, drop = True)


final_grouped_gene_df_imp_train = pd.DataFrame()

i = 0
for cell in celltype : 
	tmp_data_cell = final_grouped_gene_df[final_grouped_gene_df['CellType'] == cell]
	del tmp_data_cell['CellType']
	imp = SimpleImputer(missing_values=np.nan, strategy='mean') #python3
	tmp_data_cell_imp = imp.fit_transform(tmp_data_cell)
	tmp_data_cell_imp = pd.DataFrame(tmp_data_cell_imp)
	tmp_data_cell_imp['CellType'] = cell 
	tmp_data_cell_imp['SampleName'] = tmp_data_cell.index.tolist()
	if i == 0 :
		final_grouped_gene_df_imp_train = tmp_data_cell_imp.copy()
	else :
		final_grouped_gene_df_imp_train = pd.concat([final_grouped_gene_df_imp_train, tmp_data_cell_imp], axis = 0)
	i += 1

final_grouped_gene_df_imp_train.set_index('SampleName', inplace = True, drop = True)
final_grouped_gene_df_imp_train.columns = final_grouped_gene_df.columns

######################
####  Test  Data  ####
######################

sample_info = pd.read_csv("./testing_metadata.csv")
sample_info = sample_info.dropna(axis = 0)

filelist = pd.read_csv("testing_methyl_filelist.txt", header = None)
final_grouped_gene_df = pd.DataFrame()
for i in range(len(filelist)) : #len(filelist) 
	filename = dirname + filelist[0][i] + ".bed"
	data = pd.read_csv(filename, header = None, sep = '\t')
	del data[4]
	del data[5]
	del data[6]
	data = data.rename(columns = {0 : 'chr', 1 : 'start', 2 :'end', 3 : filelist[0][i], 7 : 'gene'})
	data_filterd = pd.merge(filtered_cell_df, data, how = 'left', on = 'gene')
	data_filterd.to_csv(resdir + filelist[0][i] + ".csv", mode = "w", index = False)
	del data_filterd['chr']
	del data_filterd['start']
	del data_filterd['end']
	grouped_gene = data_filterd.groupby('gene')
	grouped_gene_df = pd.DataFrame(grouped_gene.mean())
	if i == 0 :
		final_grouped_gene_df = grouped_gene_df.copy()
	else :
		final_grouped_gene_df = pd.merge(final_grouped_gene_df, grouped_gene_df, left_index = True, right_index = True)

final_grouped_gene_df = final_grouped_gene_df.T
tmp = final_grouped_gene_df.copy()
final_grouped_gene_df = pd.merge(final_grouped_gene_df, sample_info, left_index = True, right_on = "FileName")

del final_grouped_gene_df['OmicsType']
del final_grouped_gene_df['FileName']
final_grouped_gene_df.set_index("SampleName", inplace = True, drop = True)

train_methyl_feature = final_grouped_gene_df_imp_train.columns
test_methyl_feature = final_grouped_gene_df.columns

common_methyl_features = set(train_methyl_feature) & set(test_methyl_feature)
common_methyl_features_list = list(common_methyl_features)

final_grouped_gene_df_imp_train[common_methyl_features_list].to_csv(methyl_dir + "methyl_filtered_gene_avg_promoter_train.csv", mode = "w", index = True)
final_grouped_gene_df[common_methyl_features_list].to_csv(methyl_dir + "methyl_filtered_gene_avg_promoter_test.csv", mode = "w", index = True)

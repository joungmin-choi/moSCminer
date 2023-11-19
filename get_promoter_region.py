import pandas as pd 
import os

###########################################
### Create file list for ACC and methyl ###
###########################################

filetype_list = ['training', 'testing']
for filetype in filetype_list :
	metaData = pd.read_csv(filetype + "_metadata.csv")
	methyl = metaData[metaData['OmicsType'] == "DNA methylation"]
	acc = metaData[metaData['OmicsType'] == "DNA accessibility"]
	methyl_files = pd.DataFrame(methyl['FileName'])
	acc_files = pd.DataFrame(acc['FileName'])
	methyl_files.to_csv(filetype + "_methyl_filelist.txt", mode = "w", index = False, header = False)
	acc_files.to_csv(filetype + "_acc_filelist.txt", mode = "w", index = False, header = False)

acc_dir = "./preprocessing_acc/"
methyl_dir = "./preprocessing_methyl/"

try :
	os.mkdir(acc_dir)
except :
	print("ACC dir exist!!")

try :
	os.mkdir(methyl_dir)
except :
	print("methyl dir exist!")


'''
BED format: chr, start, end, gene_name
'''
###########################################
### Create mart_export BED file for ACC ###
###########################################

ensembl_mart_data = pd.read_csv("./mart_export.csv")
try :
	ensembl_mart_data = ensembl_mart_data.drop_duplicates(['Associated Gene Name', 'Chromosome Name', 'Transcript Start (bp)'], keep = 'first')
	ensembl_mart_data.rename(columns = {'Associated Gene Name' : 'Gene Name', 'Chromosome Name' : 'Chromosome/scaffold name', 'Transcript Start (bp)' : 'Transcript start (bp)', 'Transcript End (bp)' : 'Transcript end (bp)', 'Ensembl Gene ID' : 'Gene stable ID'}, inplace = True)
except :
	ensembl_mart_data = ensembl_mart_data.drop_duplicates(['Gene name', 'Chromosome/scaffold name', 'Transcript start (bp)'], keep = 'first')
	ensembl_mart_data.rename(columns = {'Gene name' : 'Gene Name'}, inplace = True)

ensembl_mart_data.reset_index(inplace = True, drop = True)

bed_columns_list = ['Chromosome/scaffold name', 'Transcript start (bp)', 'Transcript end (bp)', 'Gene Name']

ensembl_mart_data['Chromosome/scaffold name'] = 'chr' + ensembl_mart_data['Chromosome/scaffold name']
ensembl_mart_data_acc = ensembl_mart_data[bed_columns_list]
ensembl_mart_data_acc.to_csv(acc_dir + "mart_export_acc.bed", mode = "w", index = False, header = False, sep = '\t')

#filename = "./preprocessing_rna/ensembl_to_geneID.csv"
#gene = pd.read_csv(filename)
#gene_id = pd.DataFrame(gene["ensembl_gene_id"])
#ensembl_mart_data = pd.merge(gene_id, ensembl_mart_data, left_on = "ensembl_gene_id", right_on = "Ensembl Gene ID")
#del ensembl_mart_data['ensembl_gene_id']
#ensembl_mart_data['Chromosome Name'] = 'chr' + ensembl_mart_data['Chromosome Name']
#ensembl_mart_data_acc = ensembl_mart_data[['Chromosome Name', 'Transcript Start (bp)', 'Transcript End (bp)', 'Associated Gene Name']]
#ensembl_mart_data_acc.to_csv("./preprocessing_acc/ensembl_mart_export_common_gene.bed", mode = "w", index = False, header = False, sep = '\t')

'''
BED format: chr, start, end, gene_name
'''
##############################################
### Create mart_export BED file for methyl ###
##############################################

promoter_start = []
promoter_end = []

for i in range(len(ensembl_mart_data)) :
	if (ensembl_mart_data['Strand'][i] == -1) :
		tmp_region = ensembl_mart_data['Transcript end (bp)'][i]
		promoter_start.append(tmp_region)
		promoter_end.append(tmp_region + 2000)
	else :
		tmp_region = ensembl_mart_data['Transcript start (bp)'][i]
		if (tmp_region - 2000) < 0 :
			promoter_start.append(0)
		else :
			promoter_start.append(tmp_region - 2000)
		promoter_end.append(ensembl_mart_data['Transcript start (bp)'][i])

ensembl_mart_data_methyl = ensembl_mart_data[bed_columns_list]
ensembl_mart_data_methyl['Transcript start (bp)'] = promoter_start
ensembl_mart_data_methyl['Transcript end (bp)'] = promoter_end

ensembl_mart_data_methyl.to_csv(methyl_dir + "mart_export_methyl.bed", mode = "w", index = False, header = False, sep = '\t')









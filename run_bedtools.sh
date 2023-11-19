#!/bin/bash

acc_dir="./preprocessing_acc/"
methyl_dir="./preprocessing_methyl/"
mart_export_acc=$acc_dir"mart_export_acc.bed"
mart_export_methyl=$methyl_dir"mart_export_methyl.bed"
res_dir_acc=$acc_dir"results_find_common_region"
res_dir_methyl=$methyl_dir"results_find_common_region"

mkdir $res_dir_acc
mkdir $res_dir_methyl

dataTypeList="training testing"

for dataType in $dataTypeList
do
	acc_fileList=$dataType"_acc_filelist.txt"
	methyl_fileList=$dataType"_methyl_filelist.txt"
	data_dir="./"$dataType"_data/"

	cat $acc_fileList | while read line
	do
		echo $line
		filename=$data_dir""$line
		outputname=$line".bed"
		bedtools intersect -a $filename -b $mart_export_acc -wa -wb > $outputname
		mv $outputname $res_dir_acc
	done

	cat $methyl_fileList | while read line
	do
		echo $line
		filename=$data_dir""$line
		outputname=$line".bed"
		bedtools intersect -a $filename -b $mart_export_methyl -wa -wb > $outputname
		mv $outputname $res_dir_methyl
	done
	
done
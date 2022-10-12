#!/bin/bash
#$ -pe threads 4
#input file will be three columns delimited by a comma. {1} is the first column, {2} is the second column, {3} is the third column

#This is specific to my environment
source /mnt/home/software/anaconda/anaconda2/bin/activate ete3

#provide location of comma-separate list of files for each gene. Made with previous script. 
filename=$1
parallel --will-cite --colsep ',' /mnt/home/lgw108/.conda/envs/ete3/bin/hyphy /mnt/home/lgw108/.conda/envs/ete3/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/FitMultiModel.bf --alignment {1} --tree {2} --outfile {3} < $filename

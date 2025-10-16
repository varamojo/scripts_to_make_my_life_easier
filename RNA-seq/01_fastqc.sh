#!/bin/bash

#SBATCH -p medium
#SBATCH --job-name=fastqc
#SBATCH --output=qc_%A_%a.out
#SBATCH --error=qc_%A_%a.err
#SBATCH -c 1
#SBATCH --mem=64gb
#SBATCH --time=5:00:00
#SBATCH --array=0-24

##################################################################################################
## Script for RNAseq analysis
# Adjusted for SLURM environment, for local application, see comments below 
#
# Author: Vassiliki Varamogianni
# Created: 2025 Goettingen Germany
# Contact: vasiliki.varamogiannimamatsi@med.uni-goettingen.de
#
# Feel free to use my scripts for learning
# and non-commercial purposes. If you find it helpful, consider reaching out or giving credit!
# For commercial use, please contact the author.
#
#  DISCLAIMER: This script is provided without any warranty. Use at your own risk!
##################################################################################################

##################################################################################################
############-- EDIT HERE
##################################################################################################

ifs=(p1718sH2AZ-high-minus-ICM-minus-dox-I_hg_mrna_Papantonis_I_1/p1718sH2AZ-high-minus-ICM-minus-dox-I_Papantonis_S24_merged_L001_R1_001.fastq.gz
p1718sH2AZ-high-minus-ICM-minus-dox-I_hg_mrna_Papantonis_I_1/p1718sH2AZ-high-minus-ICM-minus-dox-I_Papantonis_S24_merged_L001_R2_001.fastq.gz
p1718sH2AZ-high-minus-ICM-minus-dox-J_hg_mrna_Papantonis_J_1/p1718sH2AZ-high-minus-ICM-minus-dox-J_Papantonis_S25_merged_L001_R1_001.fastq.gz
p1718sH2AZ-high-minus-ICM-minus-dox-J_hg_mrna_Papantonis_J_1/p1718sH2AZ-high-minus-ICM-minus-dox-J_Papantonis_S25_merged_L001_R2_001.fastq.gz
p1718sH2AZ-high-minus-ICM-plus-dox-K_hg_mrna_Papantonis_K_1/p1718sH2AZ-high-minus-ICM-plus-dox-K_Papantonis_S26_merged_L001_R1_001.fastq.gz
p1718sH2AZ-high-minus-ICM-plus-dox-K_hg_mrna_Papantonis_K_1/p1718sH2AZ-high-minus-ICM-plus-dox-K_Papantonis_S26_merged_L001_R2_001.fastq.gz
p1718sH2AZ-high-minus-ICM-plus-dox-L_hg_mrna_Papantonis_L_1/p1718sH2AZ-high-minus-ICM-plus-dox-L_Papantonis_S27_merged_L001_R1_001.fastq.gz
p1718sH2AZ-high-minus-ICM-plus-dox-L_hg_mrna_Papantonis_L_1/p1718sH2AZ-high-minus-ICM-plus-dox-L_Papantonis_S27_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-minus-ICM-minus-dox-A_hg_mrna_Papantonis_A_1/p1718sH2AZ-low-minus-ICM-minus-dox-A_Papantonis_S16_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-minus-ICM-minus-dox-A_hg_mrna_Papantonis_A_1/p1718sH2AZ-low-minus-ICM-minus-dox-A_Papantonis_S16_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-minus-ICM-minus-dox-B_hg_mrna_Papantonis_B_1/p1718sH2AZ-low-minus-ICM-minus-dox-B_Papantonis_S17_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-minus-ICM-minus-dox-B_hg_mrna_Papantonis_B_1/p1718sH2AZ-low-minus-ICM-minus-dox-B_Papantonis_S17_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-minus-ICM-plus-dox-C_hg_mrna_Papantonis_C_1/p1718sH2AZ-low-minus-ICM-plus-dox-C_Papantonis_S18_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-minus-ICM-plus-dox-C_hg_mrna_Papantonis_C_1/p1718sH2AZ-low-minus-ICM-plus-dox-C_Papantonis_S18_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-minus-ICM-plus-dox-D_hg_mrna_Papantonis_D_1/p1718sH2AZ-low-minus-ICM-plus-dox-D_Papantonis_S19_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-minus-ICM-plus-dox-D_hg_mrna_Papantonis_D_1/p1718sH2AZ-low-minus-ICM-plus-dox-D_Papantonis_S19_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-plus-ICM-minus-dox-E_hg_mrna_Papantonis_E_1/p1718sH2AZ-low-plus-ICM-minus-dox-E_Papantonis_S20_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-plus-ICM-minus-dox-E_hg_mrna_Papantonis_E_1/p1718sH2AZ-low-plus-ICM-minus-dox-E_Papantonis_S20_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-plus-ICM-minus-dox-F_hg_mrna_Papantonis_F_1/p1718sH2AZ-low-plus-ICM-minus-dox-F_Papantonis_S21_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-plus-ICM-minus-dox-F_hg_mrna_Papantonis_F_1/p1718sH2AZ-low-plus-ICM-minus-dox-F_Papantonis_S21_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-plus-ICM-plus-dox-G_hg_mrna_Papantonis_G_1/p1718sH2AZ-low-plus-ICM-plus-dox-G_Papantonis_S22_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-plus-ICM-plus-dox-G_hg_mrna_Papantonis_G_1/p1718sH2AZ-low-plus-ICM-plus-dox-G_Papantonis_S22_merged_L001_R2_001.fastq.gz
p1718sH2AZ-low-plus-ICM-plus-dox-H_hg_mrna_Papantonis_H_1/p1718sH2AZ-low-plus-ICM-plus-dox-H_Papantonis_S23_merged_L001_R1_001.fastq.gz
p1718sH2AZ-low-plus-ICM-plus-dox-H_hg_mrna_Papantonis_H_1/p1718sH2AZ-low-plus-ICM-plus-dox-H_Papantonis_S23_merged_L001_R2_001.fastq.gz)

file=${ifs[${SLURM_ARRAY_TASK_ID}]}   ## If you dont apply sbatch jobs, just replace here with a loop through your samples --> for file in ${ifs[@]};do

# export path for fastqc, md5sum and multiqc
export PATH="/mnt/vast-standard/home/varamogianni/u17810/.conda/envs/deeptools/bin/:$PATH"

# Input 
HOME="home/dir"
mkdir -p ${HOME}/RNAseq/

##################################################################################################
############-- STOP EDIT HERE !!!!!!!!!!!!!!!!!
##################################################################################################

ii=$(echo "$file" | sed -E 's/(.*_S[0-9]+).*/\1/')
echo " ${ii}"
echo "$file"
inputDir="${HOME}/RNAseq/"
outDir="${inputDir}/fastqc/"
md5Dir="inputDir/md5/"

mkdir -p ${outDir}
mkdir -p ${md5Dir}/${ii}
md5sum ${inputDir}/${file}  > ${md5Dir}/${ii}/md5sum_${file}.txt
fastqc -o ${outDir} ${inputDir}/${file}

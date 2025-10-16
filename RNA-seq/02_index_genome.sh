#!/bin/bash

#SBATCH -p medium
#SBATCH --job-name=index
#SBATCH --output=index_%A.out
#SBATCH --error=index_%A.err
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=2-00:00:00

##################################################################################################
## Script for RNAseq analysis
# Suitable for SLURM environment AND local application 
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

# load modules and export paths
module load gcc/14.2.0
module load star/2.7.11b
export PATH="/sw/rev/25.04/sapphirerapids_opx_rocky8/linux-rocky8-sapphirerapids/gcc-14.2.0/star-2.7.11b-iy3n7aikvstzl224ajgeljwxyasamuee/bin/:$PATH"

HOME="home/dir"
index="${HOME}/Mapping_and_annotation_files/star_index_2.7.11b"  ## or the path for your index 
fa="${HOME}/Mapping_and_annotation_files/references/hg38/hg38.fa"
gtf="${HOME}/Mapping_and_annotation_files/annotations/hg38/hg38.ncbiRefSeq.gtf"

##################################################################################################
############-- STOP EDIT HERE !!!!!!!!!!!!!!!!!
##################################################################################################

STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ${index} \
--genomeFastaFiles ${fa} \
--sjdbGTFfile ${gtf} \
--sjdbOverhang 100

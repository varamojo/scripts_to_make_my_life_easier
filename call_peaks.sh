#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH -c 8
#SBATCH --job-name=00_peaks
#SBATCH -o 00_peaks.%A.%a.log
#SBATCH -e 00_peaks.%A.%a.err
#SBATCH --time 03:30:00
#SBATCH --mem 10G
#SBATCH --array=0

###############################################
## Simple script for peak calling using macs2
###############################################
# Make sure you bam files are sorted and the duplicates are removed. 
# If multiple replicates, I recommend merging in the level of .bam
# See 01_bam_prepro.sh for bam handling after aligment

######################################
## START EDITING HERE!
######################################

module load macs2
module load bamtools

bamInDir="path/to/bam/"
bamInFile="file.sorted.bam"
InputDir="path/to/input"
InputFile="input.sorted.bam"
peakOut="path/to/output/peaks"
qCut="0.03"

##Call peaks using bg
macs2 callpeak -t ${bamInDir}/${bamInFile} -c ${InputDir}/${InputFile} -f BAM -g hs --bw 200 --keep-dup all -n ${peakOut}_${qCut} -q ${qCut}

##Call peaks without bg
#macs2 callpeak -t ${bamInDir}/${bamInFile} -f BAM -g hs --bw 200 --keep-dup all -n ${peakOut}_no_bg

##################################################################################################
# Author: Vassiliki Varamogianni
# Created: 2024 Goettingen Germany
# Contact: vasiliki.varamogiannimamatsi@med.uni-goettingen.de
#
# Feel free to use my scripts for learning
# and non-commercial purposes. If you find it helpful, consider reaching out or giving credit!
# For commercial use, please contact the author.
#
#  DISCLAIMER: This script is provided without any warranty. Use at your own risk!

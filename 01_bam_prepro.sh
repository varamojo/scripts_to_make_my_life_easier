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
## Simple script for .bam handling 
###############################################

######################################
## START EDITING HERE!
######################################

# Load modules
module load samtools
module load picard

#Specify parameters
Dir="path/to/working/folder/"
InBam1="bamFile1.bam"
InBam2="bamFile1.bam"
InBamMerged="bamFileMeregd.bam"

# Markduplicates
# sample 1
picard MarkDuplicates I=${Dir}/${InBam1} O=${Dir}/dedup_${InBam1} M=${Dir}/sample_metrics_${InBam1}.txt REMOVE_DUPLICATES=true
# sample 2
picard MarkDuplicates I=${Dir}/${InBam2} O=${Dir}/dedup_${InBam2} M=${Dir}/sample_metrics_${InBam2}.txt REMOVE_DUPLICATES=true


# Merge .bam files
samtools merge ${Dir}/${InBamMerged} ${Dir}/dedup_${InBam1} ${Dir}/dedup_${InBam2}

# Sort .bam files
samtools sort  ${Dir}/${InBamMerged} -o ${Dir}/Sort_${InBamMerged}

# Index .bam files
samtools index ${Dir}/Sort_${InBamMerged}




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

#!/bin/bash
# Slurm Parameters
#SBATCH -p medium
#SBATCH --job-name=puppy.all
#SBATCH -o pup_all-%A.%a.log
#SBATCH -e pup_all-%A.%a.err
#SBATCH --time 09:30:00
#SBATCH --mem 10G
#SBATCH --array=0-5

II=(04_0minLinks_NoTAD_25kb.bedpe 04_0minLinks_WithinTAD_25kb.bedpe 04_0minLinks_AcrossTADboundary_25kb.bedpe)

i=${II[${SLURM_ARRAY_TASK_ID}]}
inDir="$HOME/Links2/20240510_CoAccessibleLinks/"

mkdir -p $HOME/Links2/20240510_CoAccessibleLinks/APA_all/
outDir="$HOME/Links2/20240510_CoAccessibleLinks/APA_all/"

cool_input="$HOME/coolfiles/GSE63525_HUVEC_combined.10000kb.cool"

awk '{print $1,$2,$3,$4,$5,$6}' ${inDir}/${i} | tail -n+2 | sed 's/"//g' |  sed 's/ /\t/g' > ${inDir}/${i}.${folder}.tmp.bedpe
input_coord="${inDir}/${i}.${folder}.tmp.bedpe"

#        head ${input_coord}  --example of the format
#        1       905571  906571  1       907707  908707
#        1       910864  911864  1       917828  918828
#        1       915016  916016  1       916236  917236
#        1       925218  926218  1       931485  932485
#        1       927287  928287  1       929070  930070
#        1       929070  930070  1       930133  931133
coolpup.py ${cool_input}  ${input_coord} --outname ${outDir}/${i}.${folder}.clpy
plotpup.py ${outDir}/${i}.${folder}.clpy --output ${outDir}/${i}.${folder}.LoopGram.pdf --vmin 0.9 --vmax 1.15
# for absolute contrast
plotpup.py --input_pups ${outDir}/${i}.${folder}.clpy --output ${outDir}/${i}.${folder}.LoopGram.pdf --vmin 0.1 --vmax 1.0





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

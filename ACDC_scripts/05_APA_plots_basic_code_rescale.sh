#!/bin/bash
# Slurm Parameters
#SBATCH -p medium
#SBATCH --job-name=puppy
#SBATCH -o pup_%A.%a.log
#SBATCH -e pup_%A.%a.err
#SBATCH --time 09:30:00
#SBATCH --mem 10G
#SBATCH --array=0-9

II=(01_Arrowhead_10kb.tsv
04_Arrowhead_10kb_CoAccessibleLinksWithin.tsv
04_Arrowhead_10kb_NoCoAccessibleLink.tsv
02_Arrowhead_10kb_WithoutTSG.tsv
02_Arrowhead_10kb_WithTSG.tsv
05_Arrowhead_10kb_CoAccessibleDomainsAcross.tsv
03_Arrowhead_10kb_TSGclusterAcrossTAD.tsv
05_Arrowhead_10kb_CoAccessibleDomainsWithin.tsv
03_Arrowhead_10kb_TSGclusterWithinTAD.tsv
05_Arrowhead_10kb_NoCoAccessibleDomains.tsv
04_Arrowhead_10kb_CoAccessibleLinksAcross.tsv)

i=${II[${SLURM_ARRAY_TASK_ID}]}


dir="$HOME/Links2/20240510_TADlists_Arrowhead/"
cool_input="$HOME/coolfiles/GSE63525_HUVEC_combined.10000kb.cool"

mkdir -p ${dir}/APA_Arrowhead
outDir="${dir}/APA_Arrowhead"

#awk '{print $1,$2,$3}' ${dir}/${i} | tail -n+2 | sort -k1,1 -k2,2n | sed 's/ /\t/g' | sed 's/"//g' > ${dir}/${i}.bed
input_coord="${dir}/${i}.bed"

#    head ${input coord} --example of format
#    1       1010000 1300000
#    1       2400000 2550000
#    1       6230000 6460000
#    1       7660000 7940000
#    1       7980000 8330000

coolpup.py ${cool_input} ${input_coord} --outdir ${outDir}  --rescale --local --outname ${i}.local.rescale.clpy
plotpup.py ${outDir}/${i}.local.rescale.clpy --output  ${outDir}/${i}.local.rescale.pdf  --scale log   --vmin 0.6 --vmax 1.9

# for absolute contrast
plotpup.py ${outDir}/${i}.local.rescale.clpy --output  ${outDir}/${i}.local.rescale.pdf  --scale log   --vmin 0.6 --vmax 1.9



##################################################################################################
# Author: Vassiliki Varamogianni, MSc
# Created: 2024 Goettingen Germany
# Contact: vasiliki.varamogiannimamatsi@med.uni-goettingen.de
#
# Feel free to use my scripts for learning
# and non-commercial purposes. If you find it helpful, consider reaching out or giving credit!
# For commercial use, please contact the author.
#
#  DISCLAIMER: This script is provided without any warranty. Use at your own risk!

#!/bin/bash

# Slurm Parameters
#SBATCH -p fat
#SBATCH -c 2
#SBATCH --job-name=bwa_dove
#SBATCH -o bwa_dove_%A.%a.log
#SBATCH -e bwa_dove_%A.%a.err
#SBATCH --time 20:00:00
#SBATCH --mem 80G
#SBATCH --array=0-1

##################################################################################################
## Script for MicroC analysis, based on the pipeline provided by Dovetail
#
# Author: Vassiliki Varamogianni
# Created: 2024 Goettingen Germany
# Contact: vasiliki.varamogiannimamatsi@med.uni-goettingen.de
#
# Feel free to use my scripts for learning
# and non-commercial purposes. If you find it helpful, consider reaching out or giving credit!
# For commercial use, please contact the author.
#
#  DISCLAIMER: This script is provided without any warranty. Use at your own risk!
##################################################################################################

II=(MicroC_Huvec_Rep3_tnfa-0min MicroC_Huvec_Rep4_tnfa-30min)

#Apply samples to the parallel mode (optional)
i=${II[${SLURM_ARRAY_TASK_ID}]}

#load your tools
export PYTHONPATH=/home/.conda/envs/Dovetail/bin/:$PYTHONPATH

python -V
export PATH=$PATH://usr/users/vasiliki/bwa-0.7.16
        echo "load samtools..."
module load samtools/1.19
        echo "load bedtools..."
module load bedtools/2.29.1

# Set important directories
workdir="set/your/working/dir/"
ref_gen="$workdir/data/hg38/hg38.fa"
chr_path="$workdir/data/hg38/hg38.chrom.sizes"
fastq_dir="$workdir/fastq/samples1/"

mkdir -p $workdir/bam
output="$workdir/bam/${i}"

# Define fastq files
fastq1="${fastq_dir}/${i}_R1.fastq.gz"
fastq2="${fastq_dir}/${i}_R2.fastq.gz"

echo "${i}"

## Map to genome
bwa mem -5SP -T0 -t 24  ${ref_gen} ${fastq1} ${fastq2} -o ${output}/${i}.sam

## Recording valid ligation events
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path ${chr_path}  ${output}/${i}.sam >  ${output}/${i}.parsed.pairsam

## Sorting the pairsam file
pairtools sort --tmpdir=${output}/tmp ${output}/${i}.parsed.pairsam > ${output}/${i}.sorted.pairsam

## Removig PCR duplicates
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats ${output}/${i}.stats.txt --output ${output}/${i}.dedup.pairsam ${output}/${i}.sorted.pairsam

## Generate .pairs and bam files
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs  ${output}/${i}.mapped.pairs --output-sam ${output}/${i}.unsorted.bam  ${output}/${i}.dedup.pairsam

## Generating the final bam file
samtools sort -@16 -T ${output}/tmp/temp.bam -o ${output}/${i}.mapped.PT.bam ${output}/${i}.unsorted.bam
samtools index  ${output}/${i}.mapped.PT.bam

## Make stats file !!Download the paired get_qc.py script from the present repo!!
.conda/envs/Dovetail/bin/python3.8 get_qc.py -p ${output}/${i}.stats.txt

## Create .hic maps !!Consider downloading the relevant .jar mentioned below!!
module load r/4.0.3-java8
java -Xmx48000m -jar juicer_tools.1.9.9_jcuda.0.8.jar pre ${output}/${i}.mapped.pairs ${output}/${i}.hic  ${chr_path}  hg38

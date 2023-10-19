#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=featurecounts
#SBATCH --output=./stdout/%j_featurecounts.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/
#Version 1.5.3 on zenith 11/18/2022

gff="/home/projects-wrighton/NIH_Salmonella/all_bins_r3_r5/dRep99/dereplicated_genomes/DRAM/genes.gff"
wd="/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/"

mkdir -p featurecounts

featureCounts -T 30 -t CDS -g ID -s 2 -p -a "$gff" -o "$wd"featurecounts/counts_minid95.tsv "$wd"bowtie2_mapping/95_filtered_namesorted_bams/*.bam







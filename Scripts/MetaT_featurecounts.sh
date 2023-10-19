#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=featurecounts
#SBATCH --output=./stdout/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/projects/Salmonella_Interaction_paper/Metatranscriptomics/New_Analysis_09062022/mapping_scaffolds
#Version 1.5.3 on zenith 10/13/2022

gff="/home/projects-wrighton/NIH_Salmonella/MAGdb/all_bins_r3_r5/dRep99/dereplicated_genomes/DRAM/genes.gff"

percent_levels="90 93 95"

for perc in $percent_levels; do
  featureCounts -T 30 -t CDS -g ID -s 2 -p -a "$gff" -o ../featurecounts/counts_scaffoldmap_"$perc".tsv namesort/*_minid."$perc"_namesort.bam
done






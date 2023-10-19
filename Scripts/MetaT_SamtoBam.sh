#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=sam_to_bam
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/projects/Salmonella_Interaction_paper/Metatranscriptomics/New_Analysis_09062022/mapping_scaffolds/original_bams

sample_names="KE KG KH KT KX LK LL LM LO LP"

for name in $sample_names; do
  samtools view -@ 30 -Sb "$name".sam > "$name".bam
  rm "$name".sam
done






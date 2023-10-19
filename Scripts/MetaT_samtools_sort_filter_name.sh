#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=filter_alignment_sort
#SBATCH --output=./stdout/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/projects/Salmonella_Interaction_paper/Metatranscriptomics/New_Analysis_09062022/mapping_scaffolds

percent_levels=".97"
sample_names="KE KG KH KT KX LK LL LM LO LP"

for name in $sample_names; do
  for perc in $percent_levels; do
    reformat.sh in=original_bams/"$name".bam out="$name"_minid"$perc".bam minidfilter="$perc" primaryonly=t pairedonly=f
    samtools sort -n -@ 30 -o namesort/"$name"_minid"$perc"_namesort.bam "$name"_minid"$perc".bam
    rm "$name"_minid"$perc".bam
  done
done





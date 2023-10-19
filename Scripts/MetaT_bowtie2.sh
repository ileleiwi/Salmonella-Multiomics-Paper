#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --job-name=bowtie2_map_metaT
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/projects/Salmonella_Interaction_paper/Metatranscriptomics/New_Analysis_09062022/mapping_scaffolds

paste read1.txt read2.txt samout.txt | while read i1 i2 o1; do
  bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 40 -x ../mapping_scaffolds/bowtie_index/drepped_r3_r5_scaffolds -S "$o1" -1 "$i1" -2 "$i2"
done



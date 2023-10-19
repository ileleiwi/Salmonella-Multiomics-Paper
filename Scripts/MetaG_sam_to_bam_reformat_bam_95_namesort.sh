#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=reformat_bam_95
#SBATCH --output=./stdout/%j_reformat_bam_97.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/
#reformat.sh Written by Brian Bushnell
#Last modified January 26, 2021
#run on Zenith 11/18/22

wd="/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/"
samples="5_13_5_8 5_14_2_8 5_14_3_8 5_14_4_8 KE KL KS"

mkdir -p "$wd"bowtie2_mapping/original_bams/ "$wd"bowtie2_mapping/95_filtered_namesorted_bams/

for samp in $samples; do
  samtools view -bS -@ 30 "$wd"bowtie2_mapping/"$samp".sam > "$wd"bowtie2_mapping/original_bams/"$samp".bam
  rm "$wd"bowtie2_mapping/"$samp".sam
  reformat.sh in="$wd"bowtie2_mapping/original_bams/"$samp".bam out="$wd"bowtie2_mapping/"$samp"_filtered95.bam minidfilter=.95 primaryonly=t pairedonly=f
  samtools sort -n -@ 30 -o "$wd"bowtie2_mapping/95_filtered_namesorted_bams/"$samp"_sorted_filtered95.bam "$wd"bowtie2_mapping/"$samp"_filtered95.bam
  rm "$wd"bowtie2_mapping/"$samp"_filtered95.bam
done

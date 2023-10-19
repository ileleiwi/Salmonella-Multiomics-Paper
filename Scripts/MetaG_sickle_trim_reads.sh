#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=sickle
#SBATCH --output=./stdout/%j_sickle.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/
#sickle version 1.33
#run on Zenith 11/13/2022

mkdir -p reads/trimmed_reads reads/trimmed_reads/singles

wd="/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/"
samples="5_13_5_8 5_14_2_8 5_14_3_8 5_14_4_8 KE KL KS"
for samp in $samples; do
  sickle pe \
    -f "$wd"reads/"$samp"_R1.fastq.gz \
    -r "$wd"reads/"$samp"_R2.fastq.gz \
    -t sanger \
    -o "$wd"reads/trimmed_reads/"$samp"_trimmed_R1.fastq.gz \
    -p "$wd"reads/trimmed_reads/"$samp"_trimmed_R2.fastq.gz \
    -s "$wd"reads/trimmed_reads/singles/"$samp"_trimmed_singles.fastq \
    --g

  rm "$wd"reads/"$samp"_R1.fastq.gz "$wd"reads/"$samp"_R2.fastq.gz
done




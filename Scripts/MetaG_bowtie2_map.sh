#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=bowtie2_map
#SBATCH --output=./stdout/%j_mapbowtie2.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/
#Bowtie 2 version 2.4.5 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
#run on Zenith 11/18/22

wd="/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/"
samples="5_13_5_8 5_14_2_8 5_14_3_8 5_14_4_8 KE KL KS"

mkdir -p "$wd"bowtie2_mapping

cd "$wd"bowtie2_mapping

for samp in $samples; do
  reformat.sh -Xmx300g \
  in="$wd"reads/filtered_trimmed_reads/"$samp"/"$samp"_trimmed_R1.anqrpht.fastq.gz \
  out1="$wd"reads/filtered_trimmed_reads/"$samp"_filtered_trimmed_R1.fastq.gz \
  out2="$wd"reads/filtered_trimmed_reads/"$samp"_filtered_trimmed_R2.fastq.gz

  bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 30 \
  -x "$wd"bowtie2_db/MQHQ_scaffolds \
  -S "$wd"bowtie2_mapping/"$samp".sam \
  -1 "$wd"reads/filtered_trimmed_reads/"$samp"_filtered_trimmed_R1.fastq.gz \
  -2 "$wd"reads/filtered_trimmed_reads/"$samp"_filtered_trimmed_R2.fastq.gz
done
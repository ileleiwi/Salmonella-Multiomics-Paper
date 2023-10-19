#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=bowtie2_db
#SBATCH --output=./stdout/%j_bowtie2indexbuild.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522
#Bowtie 2 version 2.4.5 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
#run on Zenith 11/13/22

mkdir -p bowtie2_db

wd="/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/"
cd "$wd"bowtie2_db
ref_path="/home/projects-wrighton/NIH_Salmonella/all_bins_r3_r5/dRep99/dereplicated_genomes/DRAM/scaffolds.fna"
bowtie2-build "$ref_path" MQHQ_scaffolds --threads 20


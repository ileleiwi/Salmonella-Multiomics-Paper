#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=bowtie2_index_metaT
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/Salmonella/Metatranscriptomics/SalmonellaMetaTreads/New_Analysis_09062022/mapping
bowtie2-build /home/projects-wrighton/NIH_Salmonella/Salmonella/Metagenomes/MAG_database/MQ_HQ_bins/drep_out_all_final/dereplicated_genomes/mag_id/DRAM/genes.fna MQHQ_genes
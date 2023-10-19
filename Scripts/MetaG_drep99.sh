#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=99_dRep_HMP_MQHQ
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low


#wd /home/projects-wrighton/NIH_Salmonella/all_bins_r3_r5
#dRep v2.6.2
#run on Zenith 10/23/2022
dRep dereplicate dRep99 -sa 0.99 -p 30 -comp 50 -con 10 -g *.fa --genomeInfo cbajdb_r5_checkm.csv --debug

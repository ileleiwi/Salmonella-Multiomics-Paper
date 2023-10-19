#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --job-name=checkM_dram_mqhq
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/round5/metaG/pilot/all_bins
eval "$(conda shell.bash hook)"
conda activate DRAM_06132022
checkm lineage_wf -t 20  -x fa . ./checkM --tab_table
cd checkM/
checkm qa lineage.ms . -o 1 -f analyze_bins.txt --tab_table
cd ..
DRAM.py annotate -i '*.fa' -o DRAM --min_contig_size 2500 --threads 30 --checkm_quality ./checkM/analyze_bins.txt
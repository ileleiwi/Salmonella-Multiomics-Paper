#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=240:00:00
#SBATCH --mem=64gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

wd="/home/projects-wrighton/NIH_Salmonella/projects/Salmonella_Interaction_paper/Metatranscriptomics/rarify_and_remap_09242023"
samples="KE KG KH KT KX LK LL LM LO LP"
for samp in $samples; do
bbduk.sh -Xmx1G threads=50 \
ref=/opt/bbtools/bbmap/resources/phix174_ill.ref.fa.gz \
k=31 \
hdist=1 \
statscolumns=5 \
in1="$wd"/trimmed_reads/"$samp"_R1_trimmed.fastq  \
in2="$wd"/trimmed_reads/"$samp"_R2_trimmed.fastq  \
out1="$wd"/trimmed_reads/"$samp"_R1_trimmed_phix.fastq  \
out2="$wd"/trimmed_reads/"$samp"_R2_trimmed_phix.fastq  \
outm="$wd"/trimmed_reads/stats/phix.fastq.gz  \
outs="$wd"/trimmed_reads/stats/single.fastq.gz  \
refstats="$wd"/trimmed_reads/stats/phix_removed_refstats.tsv

bbmap.sh -Xmx30g threads=50 \
ref=/home/projects/Reference_Genome_Databases/JGI_contamination/ribokmers.fa.gz \
fast=f \
minid=0.90 \
local=t \
printunmappedcount=t \
in1="$wd"/trimmed_reads/"$samp"_R1_trimmed_phix.fastq \
in2="$wd"/trimmed_reads/"$samp"_R2_trimmed_phix.fastq \
outu1="$wd"/trimmed_reads/"$samp"_R1_2.5Gbps_trimmed_phix_rrna.fastq \
outu2="$wd"/trimmed_reads/"$samp"_R2_2.5Gbps_trimmed_phix_rrna.fastq \
scafstats="$wd"/trimmed_reads/stats/rRNA_2.5Gbp.tsv

rm "$wd"/trimmed_reads/"$samp"_R1_trimmed_phix.fastq "$wd"/trimmed_reads/"$samp"_R2_trimmed_phix.fastq "$wd"/trimmed_reads/"$samp"_R1_trimmed.fastq "$wd"/trimmed_reads/"$samp"_R2_trimmed.fastq
done

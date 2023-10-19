#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=rqcfilter2
#SBATCH --output=./stdout/%j_rqcfilter2.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ileleiwi@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/
#RQCFilter2 Written by Brian Bushnell
#Last modified January 26, 2021
#run on Zenith 11/13/2022
eval "$(conda shell.bash hook)"

conda activate pigz

mkdir -p reads/filtered_trimmed_reads/

wd="/home/projects-wrighton/NIH_Salmonella/Salmonella/PRJNA491522/"
samples="5_13_5_8 5_14_2_8 5_14_3_8 5_14_4_8 KE KL KS"
for samp in $samples; do
  pigz -d -p30 "$wd"reads/trimmed_reads/"$samp"_trimmed_R1.fastq.gz
  pigz -d -p30 "$wd"reads/trimmed_reads/"$samp"_trimmed_R2.fastq.gz

  rqcfilter2.sh \
    jni=t \
    threads=30 \
    in1="$wd"reads/trimmed_reads/"$samp"_trimmed_R1.fastq \
    in2="$wd"reads/trimmed_reads/"$samp"_trimmed_R2.fastq \
    path="$wd"reads/filtered_trimmed_reads/"$samp" \
    outribo="$wd"reads/filtered_trimmed_reads/rna.fq.gz \
    rna=t \
    trimfragadapter=t \
    qtrim=r \
    trimq=0 \
    maxns=1 \
    maq=10 \
    minlen=51 \
    mlf=0.33 \
    phix=t \
    removeribo=t \
    removehuman=t \
    removedog=t \
    removecat=t \
    removemouse=t \
    khist=t \
    removemicrobes=t \
    mtst=t \
    sketch=t \
    kapa=t \
    clumpify=t \
    tmpdir=null \
    barcodefilter=f \
    trimpolyg=5 \
    -Xmx300g \
    rqcfilterdata=/home/opt/RQCFilterData

done




#!/bin/bash
##SBATCH --partition=RM
##SBATCH --export=ALL
#SBATCH --partition=pfen1
#SBATCH --mem=24G
#SBATCH -t 1-0
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --signal=2
#SBATCH --job-name=rheMac8
#SBATCH --error=logs/rheMac8-index_genome_%A_out.txt
#SBATCH --output=logs/rheMac8-index_genome_%A_out.txt

#################################
## make a bowtie2 rheMac8 index 
SCRATCH=$HOME
PROJDIR=$SCRATCH/projects/macaque_snATAC-seq
GENOMEDIR=$SCRATCH/resources/genomes

## get the UCSC rheMac8 genome
mkdir -p $GENOMEDIR/rheMac8
cd $GENOMEDIR/rheMac8

faidx $GENOMEDIR/rheMac8/rheMac8.fa -i chromsizes > rheMac8_sizes.chrom

source activate snap
ALIGNER_PATH=$(dirname -- "$(which bowtie2)")/

snaptools index-genome \
--input-fasta=$GENOMEDIR/rheMac8/rheMac8.fa \
--output-prefix=rheMac8 \
--aligner=bowtie2 \
--path-to-aligner=$ALIGNER_PATH \
--num-threads=8

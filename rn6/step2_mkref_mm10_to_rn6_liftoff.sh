#!/bin/bash
#SBATCH --partition=RM
#SBATCH -t 1-0
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --job-name=rn6
#SBATCH --error=logs/rn6_genome_%A_out.txt
#SBATCH --output=logs/rn6_genome_%A_out.txt

##################################################
## using GENCODE human gene annotation to guide the rn6 gene annotation
## using liftoff rather than liftover 
## see: https://github.com/agshumate/Liftoff

GENOME=rn6
GENOMEDIR=/home/bnphan/resources/genomes
TMPDIR=/scratch/bnphan/$GENOME
TARGET_GENOME=${GENOMEDIR}/$GENOME/$GENOME.fa
REFERENCE_GENOME=/home/bnphan/resources/genomes/mm10/mm10.fa
REFERENCE_ANNOT=/home/bnphan/resources/genomes/mm10/mm10.ncbiRefSeq.gtf

mkdir -p $GENOMEDIR/$GENOME $TMPDIR
cd $GENOMEDIR/$GENOME

# ##  download the refseq mm10 genome sequence 
if [ ! -f $TARGET_GENOME ]; then
cd $GENOMEDIR/$GENOME
wget https://hgdownload.soe.ucsc.edu/goldenPath/$GENOME/bigZips/$GENOME.fa.gz
gunzip $GENOME.fa.gz
fi

# ##  download the refseq mm10 genome sequence 
if [ ! -f $REFERENCE_GENOME ]; then
cd /home/bnphan/resources/genomes/mm10
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
gunzip mm10.ncbiRefSeq.gtf.gz
fi

################################
## liftOff the RefSeq human annotations to the rn6 genome

# generate a liftoff gene annotation from mm10 to rn6
liftoff -infer_genes -dir $TMPDIR -g $REFERENCE_ANNOT \
-o $GENOMEDIR/$GENOME/${GENOME}_liftoff_mm10_RefSeq.gff3 \
$TARGET_GENOME $REFERENCE_GENOME

# convert gff to gtf, and filter on biotypes
gffread rn6_liftoff_mm10_RefSeq.gff3 -FT -o rn6_liftoff_mm10_RefSeq.gtf

## liftover mm10 blacklist to rn6
liftOver mm10-blacklist.v2.bed \
/home/bnphan/resources/liftOver_chainz/mm10ToRn6.over.chain.gz \
rn6_liftOver_mm10-blacklist.v2.bed unlifted.bed

rm mm10-blacklist.v2.bed




#!/bin/bash
#SBATCH --partition=RM
#SBATCH -t 1-0
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --job-name=rheMac8
#SBATCH --error=logs/rheMac8_genome_%A_out.txt
#SBATCH --output=logs/rheMac8_genome_%A_out.txt

##################################################
## make a cellranger-atac reference genome for rheMac8
##
## using GENCODE human gene annotation to guide the rheMac8 gene annotation
## using liftoff rather than liftover 
## see: https://github.com/agshumate/Liftoff
## 
## use 

CELLRANGERDIR=$SCRATCH/src/cellranger-atac-1.2.0
GENOMEDIR=~/resources/genomes

########################
# get the human genome
mkdir -p $GENOMEDIR/GRCh38.p13
cd $GENOMEDIR/GRCh38.p13

## download the hg38 blacklist from ENCODE
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz

## download the hg38 to rheMac8 chain
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToRheMac8.over.chain.gz

################################
## get the UCSC rheMac8 genome
mkdir -p $GENOMEDIR/rheMac8
cd $GENOMEDIR/rheMac8

# download the rheMac8 genome sequences 
wget http://hgdownload.cse.ucsc.edu/goldenPath/rheMac8/bigZips/rheMac8.fa.gz
gunzip rheMac8.fa.gz
faidx rheMac8.fa -i bed > rheMac8_sizes.bed

# download the rheMac8 ensemble gene annotation
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/genes/rheMac8.ensGene.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/genes/rheMac8.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/genes/rheMac8.refGene.gtf.gz
gunzip rheMac8.ensGene.gtf.gz
gunzip rheMac8.ncbiRefSeq.gtf.gz
gunzip rheMac8.refGene.gtf.gz

## liftover hg38 blacklist regions to rheMac8
liftOver \
$GENOMEDIR/GRCh38.p13/hg38-blacklist.v2.bed \
$GENOMEDIR/GRCh38.p13/hg38ToRheMac8.over.chain.gz \
tmp.bed unmapped_hg38_blacklist.txt
cut -f 1,2,3 tmp.bed > rheMac8_liftOver_hg38-blacklist.v2.bed
rm tmp.bed

# generate a liftoff gene annotation from GRCh38 RefSeq NCBI to rheMac8
liftoff \
-g $GENOMEDIR/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff \
-o $GENOMEDIR/rheMac8/rheMac8_liftoff_GRCh38.p13_RefSeq.gff3 \
$GENOMEDIR/rheMac8/rheMac8.fa \
$GENOMEDIR/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna

# convert gff to gtf
gffread $GENOMEDIR/rheMac8/rheMac8_liftoff_GRCh38.p13_RefSeq.gff3 \
-FT -o $GENOMEDIR/rheMac8/rheMac8_liftoff_GRCh38.p13_RefSeq.gtf



#################################################
## create the reference genome for cellranger-atac
cd ${GENOMEDIR}

# make ref genome w/ rheMac8 liftoff annotations
REFGENOME1='refdata-cellranger-atac-rheMac8_liftoff_GRCh38.p13-1.2.0'
rm -r $REFGENOME1
${CELLRANGERDIR}/cellranger-atac mkref $REFGENOME1 \
	--config $PROJDIR/genomes/rheMac8/rheMac8_liftoff_GRCh38.p13_config.txt

# make ref genome w/ rheMac8 ensGene annotations
REFGENOME2='refdata-cellranger-atac-rheMac8-1.2.0'
rm -r $REFGENOME2
${CELLRANGERDIR}/cellranger-atac mkref $REFGENOME2 \
	--config $PROJDIR/genomes/rheMac8/rheMac8_config.txt







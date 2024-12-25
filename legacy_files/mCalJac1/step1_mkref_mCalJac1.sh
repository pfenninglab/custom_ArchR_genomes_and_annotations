#!/bin/bash
#SBATCH --partition=RM
#SBATCH -t 1-0
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --job-name=mCalJac1
#SBATCH --error=logs/mCalJac1_genome_%A_out.txt
#SBATCH --output=logs/mCalJac1_genome_%A_out.txt

##################################################
## make a cellranger-atac reference genome for mCalJac1
##
## using GENCODE human gene annotation to guide the mCalJac1 gene annotation
## using liftoff rather than liftover 
## see: https://github.com/agshumate/Liftoff

GENOMEDIR=/home/bnphan/resources/genomes
TARGET_GENOME1=${GENOMEDIR}/mCalJac1/GCA_011078405.1_mCalJac1.mat_genomic.fna
TARGET_GENOME2=${GENOMEDIR}/mCalJac1/GCA_011078405.1_mCalJac1.mat_genomic.UCSCchr.fna
REFERENCE_GENOME=/home/bnphan/resources/genomes/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna
REFERENCE_ANNOT=/home/bnphan/resources/genomes/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff

# ##  download the refseq GRCh38 genome sequence 
if [ ! -f $REFERENCE_GENOME ]; then
cd /home/bnphan/resources/genomes/GRCh38.p13
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
fi

################################
## liftOff the RefSeq human annotations to the mCalJac1 genome
mkdir -p $GENOMEDIR/mCalJac1
cd $GENOMEDIR/mCalJac1

#################################################################
## change to UCSC chromosome naming convention in mCalJac1 genome
awk -v "OFS=\t" 'NR>31 {if ( $5 ~ /CM/ ) {$3 = "chr"$3} else {$3=$5}; print $5, $3}' GCA_011078405.1_mCalJac1.mat_assembly_report.txt > GCA_011078405.1_mCalJac1.ucscChr.txt
sed -e $(awk '{print "s/"$1"/"$2"/"}' GCA_011078405.1_mCalJac1.ucscChr.txt |tr '\n' ';' ) $TARGET_GENOME1 > $TARGET_GENOME2

# generate a liftoff gene annotation from GRCh38 to mCalJac1
liftoff -g $REFERENCE_ANNOT \
-o $GENOMEDIR/mCalJac1/mCalJac1_liftoff_GRCh38.p13_RefSeq.gff3 \
$TARGET_GENOME2 $REFERENCE_GENOME

# convert gff to gtf, and filter on biotypes
gffread mCalJac1_liftoff_GRCh38.p13_RefSeq.gff3 -FT -o mCalJac1_liftoff_GRCh38.p13_RefSeq.gtf
cellranger mkgtf $GENOMEDIR/mCalJac1/mCalJac1_liftoff_GRCh38.p13_RefSeq.gtf \
$GENOMEDIR/mCalJac1/mCalJac1_liftoff_GRCh38.p13_RefSeq.filtered.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lincRNA \
--attribute=gene_biotype:antisense \
--attribute=gene_biotype:IG_LV_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:IG_V_pseudogene \
--attribute=gene_biotype:IG_D_gene \
--attribute=gene_biotype:IG_J_gene \
--attribute=gene_biotype:IG_J_pseudogene \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_C_pseudogene \
--attribute=gene_biotype:TR_V_gene \
--attribute=gene_biotype:TR_V_pseudogene \
--attribute=gene_biotype:TR_D_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_J_pseudogene \
--attribute=gene_biotype:TR_C_gene
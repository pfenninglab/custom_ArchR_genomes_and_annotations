#!/bin/bash
#SBATCH --partition=RM
#SBATCH -t 1-0
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --signal=2
#SBATCH --job-name=gffcompare
#SBATCH --dependency=afterok:10624939
#SBATCH --error=logs/rheMac10_compare_%A_out.txt
#SBATCH --output=logs/rheMac10_compare_%A_out.txt


##################################################
## make a cellranger-atac reference genome for rheMac10
##
## using GENCODE human gene annotation to guide the rheMac10 gene annotation
## using liftoff rather than liftover 
## see: https://github.com/agshumate/Liftoff
## 
## use 

PROJDIR=$SCRATCH/projects/macaque_snATAC-seq
GENOMEDIR=$SCRATCH/resources/genomes
PREMRNA=/projects/pfenninggroup/machineLearningForComputationalBiology/SingleCell/genomes/rhemac10.premrna.gtf

mkdir -p $PROJDIR/genomes/gffcompare
cd $PROJDIR/genomes/gffcompare

###############################################
## compare gene annotations across official 

## rheMac10 annotations vs. GRCh38.p13 liftoff 
gffcompare -r $GENOMEDIR/rheMac10/rheMac10.refGene.gtf \
-s $GENOMEDIR/rheMac10/rheMac10.fa -V -T \
-o rheMac10_refGene_vs_liftoff_GRCh38 \
$GENOMEDIR/rheMac10/rheMac10_liftoff_GRCh38.p13_RefSeq.gff3 

gffcompare -r $GENOMEDIR/rheMac10/rheMac10.ensGene.gtf \
-s $GENOMEDIR/rheMac10/rheMac10.fa -V -T \
-o rheMac10_ensGene_vs_liftoff_GRCh38 \
$GENOMEDIR/rheMac10/rheMac10_liftoff_GRCh38.p13_RefSeq.gff3 

gffcompare -r $GENOMEDIR/rheMac10/rheMac10.ncbiRefSeq.gtf \
-s $GENOMEDIR/rheMac10/rheMac10.fa -V -T \
-o rheMac10_ncbiRefSeq_vs_liftoff_GRCh38 \
$GENOMEDIR/rheMac10/rheMac10_liftoff_GRCh38.p13_RefSeq.gff3 


## using the rheMac10_liftoff_GRCh38.p13 annotation as reference
gffcompare -r $GENOMEDIR/rheMac10/rheMac10_liftoff_GRCh38.p13_RefSeq.gff3 \
-s $GENOMEDIR/rheMac10/rheMac10.fa -V -T \
-o rheMac10_liftoff_GRCh38_vs_refGeneEnsGeneNcbiRefseq \
$GENOMEDIR/rheMac10/rheMac10.refGene.gtf \
$GENOMEDIR/rheMac10/rheMac10.ensGene.gtf \
$GENOMEDIR/rheMac10/rheMac10.ncbiRefSeq.gtf

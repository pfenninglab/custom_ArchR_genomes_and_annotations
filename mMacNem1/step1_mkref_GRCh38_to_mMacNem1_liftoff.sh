#!/bin/bash
#SBATCH --partition=RM
#SBATCH -t 1-0
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --job-name=mMacNem1
#SBATCH --error=logs/mMacNem1_genome_%A_out.txt
#SBATCH --output=logs/mMacNem1_genome_%A_out.txt

##################################################
## using GENCODE human gene annotation to guide the mMacNem1 gene annotation
## using liftoff rather than liftover 
## see: https://github.com/agshumate/Liftoff

GENOME=mMacNem1
GENOMEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/resources/genomes
LIFTOVER=/projects/pfenninggroup/machineLearningForComputationalBiology/resources/liftOver_chainz
TMPDIR=/scratch/bnphan/$GENOME
TARGET_GENOME=${GENOMEDIR}/$GENOME/mMacNem1.hap1.cur.20240610.fasta
REFERENCE_GENOME=/projects/pfenninggroup/machineLearningForComputationalBiology/resources/genomes/GRCh38.p13/GRCh38.primary_assembly.genome.fa
REFERENCE_ANNOT=/home/bnphan/src/cellranger-6.0.0/genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf

mkdir -p $GENOMEDIR/$GENOME $TMPDIR
cd $GENOMEDIR/$GENOME

# ##  download the refseq GRCh38 genome sequence 
if [ ! -f $TARGET_GENOME ]; then
cd $GENOMEDIR/$GENOME
wget https://s3.amazonaws.com/genomeark/species/Macaca_nemestrina/mMacNem1/assembly_curated/mMacNem1.hap1.cur.20240610.fasta.gz
gunzip $GENOME.gz
fi

################################################################################################
## liftOff the CellRanger-version of GENCODE human gene annotations to the mMacNem1 genome
# generate a liftoff gene annotation from GRCh38 to mMacNem1
mamba activate liftoff

liftoff -dir $TMPDIR -g $REFERENCE_ANNOT \
-o $GENOMEDIR/$GENOME/${GENOME}_liftoff_GRCh38.gff3 \
$TARGET_GENOME $REFERENCE_GENOME

conda deactivate

# convert gff to gtf, and filter on biotypes
gffread mMacNem1_liftoff_GRCh38.gff3 -FT -o mMacNem1_liftoff_GRCh38.gtf

gzip *.gtf *.gff3




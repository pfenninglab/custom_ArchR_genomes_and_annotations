#!/bin/bash
#SBATCH --partition=RM
#SBATCH -t 1-0
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --job-name=macFas5
#SBATCH --error=logs/macFas5_genome_%A_out.txt
#SBATCH --output=logs/macFas5_genome_%A_out.txt

##################################################
## using GENCODE human gene annotation to guide the macFas5 gene annotation
## using liftoff rather than liftover 
## see: https://github.com/agshumate/Liftoff

GENOME=macFas5
GENOMEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/resources/genomes
LIFTOVER=/projects/pfenninggroup/machineLearningForComputationalBiology/resources/liftOver_chainz
TMPDIR=/scratch/bnphan/$GENOME
TARGET_GENOME=${GENOMEDIR}/$GENOME/$GENOME.fa
REFERENCE_GENOME=/projects/pfenninggroup/machineLearningForComputationalBiology/resources/genomes/GRCh38.p13/GRCh38.primary_assembly.genome.fa
REFERENCE_ANNOT=/home/bnphan/src/cellranger-6.0.0/genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf

mkdir -p $GENOMEDIR/$GENOME $TMPDIR
cd $GENOMEDIR/$GENOME

# ##  download the refseq GRCh38 genome sequence 
if [ ! -f $TARGET_GENOME ]; then
cd $GENOMEDIR/$GENOME
wget https://hgdownload.soe.ucsc.edu/goldenPath/$GENOME/bigZips/$GENOME.fa.gz
gunzip $GENOME.fa.gz
fi

################################################################################################
## liftOff the CellRanger-version of GENCODE human gene annotations to the macFas5 genome

conda activate liftoff

# generate a liftoff gene annotation from GRCh38 to macFas5
liftoff -dir $TMPDIR -g $REFERENCE_ANNOT \
-o $GENOMEDIR/$GENOME/${GENOME}_liftoff_GRCh38.gff3 \
$TARGET_GENOME $REFERENCE_GENOME

conda deactivate

# convert gff to gtf, and filter on biotypes
gffread macFas5_liftoff_GRCh38.gff3 -FT -o macFas5_liftoff_GRCh38.gtf


# ##  download the refseq GRCh38 genome sequence 
if [ ! -f $LIFTOVER/hg38ToMacFas5.over.chain.gz ]; then
cd $LIFTOVER
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMacFas5.over.chain.gz
cd $GENOMEDIR/$GENOME
fi


############################################
## liftover GRCh38 blacklist to macFas5 ####
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz

liftOver GRCh38-blacklist.v2.bed \
/projects/pfenninggroup/machineLearningForComputationalBiology/resources/liftOver_chainz/hg38ToMacFas5.over.chain.gz \
macFas5_liftOver_GRCh38-blacklist.v2.bed unlifted.bed

rm GRCh38-blacklist.v2.bed






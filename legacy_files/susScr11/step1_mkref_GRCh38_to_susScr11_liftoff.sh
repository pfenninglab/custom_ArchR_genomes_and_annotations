#!/bin/bash
#SBATCH --partition=RM
#SBATCH -t 1-0
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --job-name=susScr11
#SBATCH --error=logs/susScr11_genome_%A_out.txt
#SBATCH --output=logs/susScr11_genome_%A_out.txt

##################################################
## using GENCODE human gene annotation to guide the susScr11 gene annotation
## using liftoff rather than liftover 
## see: https://github.com/agshumate/Liftoff

GENOME=susScr11
GENOMEDIR=/home/bnphan/resources/genomes
TMPDIR=/scratch/bnphan/$GENOME
TARGET_GENOME=${GENOMEDIR}/$GENOME/$GENOME.fa
REFERENCE_GENOME=/home/bnphan/resources/genomes/GRCh38.p13/GRCh38.primary_assembly.genome.fa
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
## liftOff the CellRanger-version of GENCODE human gene annotations to the susScr11 genome

conda activate liftoff

# generate a liftoff gene annotation from GRCh38 to susScr11
liftoff -dir $TMPDIR -g $REFERENCE_ANNOT \
-o $GENOMEDIR/$GENOME/${GENOME}_liftoff_GRCh38.gff3 \
$TARGET_GENOME $REFERENCE_GENOME

conda deactivate

# convert gff to gtf, and filter on biotypes
gffread susScr11_liftoff_GRCh38.gff3 -FT -o susScr11_liftoff_GRCh38.gtf

## liftover GRCh38 blacklist to susScr11
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz

liftOver GRCh38-blacklist.v2.bed \
/home/bnphan/resources/liftOver_chainz/hg38ToSusScr11.over.chain.gz \
susScr11_liftOver_GRCh38-blacklist.v2.bed unlifted.bed

rm GRCh38-blacklist.v2.bed






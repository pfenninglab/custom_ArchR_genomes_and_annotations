GENOME_FASTA=/home/bnphan/resources/genomes/mCalJac1/GCA_011078405.1_mCalJac1.mat_genomic.UCSCchr.fna
GENE_ANNOT=/home/bnphan/resources/genomes/mCalJac1/mCalJac1_liftoff_GRCh38.p13_RefSeq.filtered.gtf
GENOMEDIR=/home/bnphan/resources/genomes/mCalJac1
mkdir -p $GENOMEDIR; cd $GENOMEDIR

## Make the genome index for bowtie2
# conda activate snap
# bowtie2-build --large-index --threads 8 $GENOME_FASTA mCalJac1
# hisat2-build -p 8 $GENOME_FASTA mCalJac1

# ####################################
# ## 2) Make the genome index for STAR
# GENE_ANNOT3=$GENOMEDIR/mCalJac1.Human_GenBankNames_modifiedPlus.filtered.gff

# awk  -v OFS='\t' '{ if ($2 == "HAVANA" || $2 == "ENSEMBL" ) { print } }' $GENE_ANNOT > $GENE_ANNOT2

# awk  -v OFS='\t' '{ if ($2 == "HAVANA" || $2 == "ENSEMBL") { 
# 	if ( $9 ~ /gene_type=protein_coding/ || $9 ~ /gene_type=lincRNA/ || $9 ~ /gene_type=antisense/ )
# 	 { print } } }' $GENE_ANNOT > $GENE_ANNOT3

STAR --runThreadN 14 \
--runMode genomeGenerate \
--genomeDir $GENOMEDIR \
--genomeFastaFiles $GENOME_FASTA \
--sjdbGTFfile $GENE_ANNOT \
--sjdbOverhang 100

mkdir -p $GENOMEDIR/mCalJac1_cellrangerRef
cd $GENOMEDIR/mCalJac1_cellrangerRef

cellranger mkref \
--genome=mCalJac1 \
--fasta=$GENOME_FASTA \
--genes=$GENE_ANNOT
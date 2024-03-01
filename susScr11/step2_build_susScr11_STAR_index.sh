GENOMEDIR=/home/bnphan/resources/genomes/susScr11
GENOME_FASTA=${GENOMEDIR}/susScr11.fa
GENE_ANNOT=${GENOMEDIR}/susScr11_liftoff_GRCh38.gtf
mkdir -p $GENOMEDIR; cd $GENOMEDIR

####################################
## 1) Make the genome index for STAR
STAR --runThreadN 14 \
--runMode genomeGenerate \
--genomeDir $GENOMEDIR \
--genomeFastaFiles $GENOME_FASTA \
--sjdbGTFfile $GENE_ANNOT \
--sjdbOverhang 100

GENOMEDIR=/home/bnphan/resources/genomes/rheMac10
GENOME_FASTA=${GENOMEDIR}/rheMac10.fa
GENE_ANNOT=${GENOMEDIR}/rheMac10_liftoff_GRCh38.p13_RefSeq.gff3
mkdir -p $GENOMEDIR; cd $GENOMEDIR

if [ ! -f $GENOME_FASTA ]; then
echo 'Getting rheMac10 genome fasta.'
wget http://hgdownload.cse.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz
gunzip rheMac10.fa.gz
fi 

####################################
## 1) Make the genome index for STAR
GENE_ANNOT2=${GENOMEDIR}/rheMac10_liftoff_GRCh38.p13_RefSeq.gtf
if [[ ! -f $GENE_ANNOT2 ]]; then gffread $GENE_ANNOT -E -F -T -o $GENE_ANNOT2; fi
head -3000 $GENE_ANNOT2 > tmp.txt

STAR --runThreadN 14 \
--runMode genomeGenerate \
--genomeDir $GENOMEDIR \
--genomeFastaFiles $GENOME_FASTA \
--sjdbGTFfile $GENE_ANNOT2 \
--sjdbOverhang 100
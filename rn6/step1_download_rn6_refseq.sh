
PROJDIR=/home/bnphan/resources/genomes/rn6
cd $PROJDIR

## download from https://ftp.cngb.org/pub/CNSA/data3/CNP0001471/
wget -r --no-parent -nH --cut-dirs=10 --tries=0 \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/895/GCA_000001895.4_Rnor_6.0/

wget https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/genes/rn6.ncbiRefSeq.gtf.gz

awk -F"," 'FNR==NR + 35{var[$5]=$1;next;}{print var[$1]FS$2}'  GCA_000001895.4_Rnor_6.0_assembly_report.txt file2

## download the mouse blacklist
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
gunzip mm10-blacklist.v2.bed.gz

liftOver mm10-blacklist.v2.bed \
/home/bnphan/resources/liftOver_chainz/mm10ToRn6.over.chain.gz \
rn6_liftOver_mm10-blacklist.v2.bed unlifted.bed

rm mm10-blacklist.v2.bed

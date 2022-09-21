ss <- function(x, pattern, slot = 1, ...) {
  sapply(strsplit(x = x, split = pattern, ...), '[', slot)
}
if(FALSE){
  # install the rheMac8 references
  BiocManager::install("BSgenome.Mmulatta.UCSC.rheMac8")
  BiocManager::install("org.Mmu.eg.db")
} 

###########################################
## load in all the required packages ######
suppressPackageStartupMessages(library(ArchR))
library(BSgenome.Mmulatta.UCSC.rheMac8)
library(GenomicFeatures)
library(AnnotationDbi)

GENOMEDIR='/home/bnphan/resources/genomes/rheMac8'

###########################################################
###### create the genome annotation files for ArchR #######
genomeAnnotation = createGenomeAnnotation(genome = BSgenome.Mmulatta.UCSC.rheMac8, 
                                          filter = TRUE, filterChr = c("chrM"))
chromSizes = genomeAnnotation$chromSizes
genome(chromSizes) <- "rheMac8"
start(chromSizes) = start(chromSizes)+1e6
end(chromSizes) = end(chromSizes)-1e6

blacklist = import(file.path(GENOMEDIR,'rheMac8_liftOver_hg38-blacklist.v2.bed'))
blacklist = blacklist[!grepl('NW',seqnames(blacklist))]
blacklist <- sort(sortSeqlevels(blacklist), ignore.strand = TRUE)
seqlevels(blacklist) <- seqlevels(chromSizes)
seqlengths(blacklist) = seqlengths(chromSizes)
genome(blacklist) <- "rheMac8"

blacklist = trim(blacklist)
genomeAnnotation$blacklist = blacklist


#################################################################
###### load in all the genomes, annotations, for rheMac8 #######
if(FALSE){
  # load in the gene annotation and save to sqlite 
  txdb = makeTxDbFromGFF(
    file.path(GENOMEDIR,'rheMac8_liftoff_GRCh38.p13_RefSeq.gff3'), 
    organism = 'Macaca mulatta', dbxrefTag = 'Dbxref') # for Macaca mulatta
  seqlevels(txdb) <- seqlevels(chromSizes)
  saveDb(txdb, file.path(GENOMEDIR,'rheMac8_liftoff_GRCh38.p13_RefSeq.sqlite'))
  
  # format genes
  library(org.Hs.eg.db) # the liftoff gene symbols from hg38 by ENSEMBL ID
  genes <- GenomicFeatures::genes(txdb)
  genes = genes[ !duplicated(genes)]
  seqlevels(genes) <- seqlevels(genes)[!grepl('NW',seqlevels(genes))]
  mcols(genes)$symbol <- mcols(genes)$gene_id
  genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
  seqlevels(genes, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(genes) = end(chromSizes)
  genome(genes) <- "rheMac8"
  genes = trim(genes)
  
  # format exons
  exons <- unlist(GenomicFeatures::exonsBy(txdb, by = "tx"))
  seqlevels(exons) <- seqlevels(exons)[!grepl('NW',seqlevels(exons))]
  exons$tx_id <- names(exons)
  mcols(exons)$symbol <- select(txdb, keys = paste0(mcols(exons)$tx_id), 
                                column = "GENEID", keytype = "TXID")[, "GENEID"]
  names(exons) <- NULL
  mcols(exons)$exon_id <- NULL
  mcols(exons)$exon_name <- NULL
  mcols(exons)$exon_rank <- NULL
  mcols(exons)$tx_id <- NULL
  exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)
  exons = exons[!is.na(exons$symbol) & !duplicated(exons)]
  seqlevels(exons, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(exons) = end(chromSizes)
  genome(exons) <- "rheMac8"
  exons = trim(exons)
  
  # TSS
  TSS <- resize(genes, 1, "start")
  TSS <- sort(sortSeqlevels(TSS), ignore.strand = TRUE)
  seqlevels(TSS, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(TSS) = end(chromSizes)
  genome(TSS) <- "rheMac8"
  TSS = trim(TSS)
  
  # save the files
  save(genes, exons, TSS, file = 
  file.path(GENOMEDIR,'rheMac8_liftoff_GRCh38.p13_RefSeq_genes_exons_TSS.rda'))
} else {
  txdb = loadDb(file=file.path(GENOMEDIR,'rheMac8_liftoff_GRCh38.p13_RefSeq.sqlite'))
  load(file.path(GENOMEDIR,'rheMac8_liftoff_GRCh38.p13_RefSeq_genes_exons_TSS.rda'))
}

# make ArchR gene annotation 
geneAnnotation <- createGeneAnnotation(genes = genes, exons = exons, TSS = TSS)

save(genomeAnnotation, geneAnnotation, file = 
       file.path(GENOMEDIR,'rheMac8_liftoff_GRCh38.p13_ArchR_annotations.rda'))


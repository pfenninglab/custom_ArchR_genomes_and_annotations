ss <- function(x, pattern, slot = 1, ...) {
  sapply(strsplit(x = x, split = pattern, ...), '[', slot)
}

###########################################
## load in all the required packages ######
suppressPackageStartupMessages(library(ArchR))
library(BSgenome.mCalJac1.mat)
library(GenomicFeatures)
library(AnnotationDbi)

GENOMEDIR='/home/bnphan/resources/genomes/mCalJac1'

###########################################################
###### create the genome annotation files for ArchR #######
genomeAnnotation = createGenomeAnnotation(
  genome = BSgenome.mCalJac1.mat, filter = TRUE, 
  filterChr = grep('JAA', seqnames(BSgenome.mCalJac1.mat), value = T))

chromSizes = genomeAnnotation$chromSizes
genome(chromSizes) <- "mCalJac1"

# to avoid ends of chromosomes tiled regions
restrict = chromSizes
start(restrict) = start(restrict)+1e6
end(restrict) = end(restrict)-1e6

#################################################################
###### load in all the genomes, annotations, for mCalJac1 #######
if(FALSE){
  # load in the gene annotation and save to sqlite 
  txdb = makeTxDbFromGFF(
    file.path(GENOMEDIR,'mCalJac1_liftoff_GRCh38.p13_RefSeq.gff3'), 
    organism = 'Callithrix jacchus', dbxrefTag = 'Dbxref') # for Callithrix jacchus
  seqlevels(txdb) <- seqlevels(chromSizes)
  saveDb(txdb, file.path(GENOMEDIR,'mCalJac1_liftoff_GRCh38.p13_RefSeq.sqlite'))
}else {
  txdb = loadDb(file.path(GENOMEDIR,'mCalJac1_liftoff_GRCh38.p13_RefSeq.sqlite'))
}  

if(FALSE){
  # format genes
  library(org.Hs.eg.db) # the liftoff gene symbols from hg38 by ENSEMBL ID
  genes <- GenomicFeatures::genes(txdb)
  genes = genes[ !duplicated(genes)]
  genes <- dropSeqlevels(genes, grep('NW',seqlevels(genes), value = T),
                         pruning.mode="coarse")
  mcols(genes)$symbol <- mcols(genes)$gene_id
  genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
  seqlevels(genes, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(genes) = end(chromSizes)
  genome(genes) <- "mCalJac1"
  genes = subsetByOverlaps(genes, restrict)
  
  # format exons
  exons <- unlist(GenomicFeatures::exonsBy(txdb, by = "tx"))
  exons <- dropSeqlevels(exons, grep('NW',seqlevels(exons), value = T),
                         pruning.mode="coarse")
  exons$tx_id <- names(exons)
  mcols(exons)$symbol <- select(txdb, keys = paste0(mcols(exons)$tx_id), 
                                column = "GENEID", keytype = "TXID")[, "GENEID"]
  names(exons) <- NULL
  mcols(exons)$exon_id <- NULL
  mcols(exons)$exon_name <- NULL
  mcols(exons)$exon_rank <- NULL
  mcols(exons)$tx_id <- NULL
  exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)
  exons = exons[!is.na(exons$symbol) & !duplicated(exons) & 
                  exons$symbol %in% genes$symbol]
  seqlevels(exons, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(exons) = end(chromSizes)
  genome(exons) <- "mCalJac1"
  genes = genes[genes$symbol %in% unique(exons$symbol)]
  # exons = subsetByOverlaps(exons, restrict)
  
  # TSS
  TSS <- resize(genes, 1, "start")
  TSS <- sort(sortSeqlevels(TSS), ignore.strand = TRUE)
  TSS <- dropSeqlevels(TSS, grep('NW',seqlevels(TSS), value = T),
                       pruning.mode="coarse")
  seqlengths(TSS) = end(chromSizes)
  genome(TSS) <- "mCalJac1"
  # TSS = subsetByOverlaps(TSS, restrict)
  
  # save the files
  save(genes, exons, TSS, file = 
         file.path(GENOMEDIR,'mCalJac1_liftoff_GRCh38.p13_RefSeq_genes_exons_TSS.rda'))
} else {
  load(file.path(GENOMEDIR,'mCalJac1_liftoff_GRCh38.p13_RefSeq_genes_exons_TSS.rda'))
}

# make ArchR gene annotation 
geneAnnotation <- createGeneAnnotation(genes = genes, exons = exons, TSS = TSS)
save(genomeAnnotation, geneAnnotation, file = 
       file.path(GENOMEDIR,'mCalJac1.mat_liftoff_GRCh38.p13_ArchR_annotations.rda'))

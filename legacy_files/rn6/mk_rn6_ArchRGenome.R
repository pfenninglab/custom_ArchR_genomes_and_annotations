ss <- function(x, pattern, slot = 1, ..) {
  sapply(strsplit(x = x, split = pattern, ..), '[', slot)
}
if(FALSE){
  # install the rn6 references
  BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn6")
  BiocManager::install("org.Rn.eg.db")
} 

###########################################
## load in all the required packages ######
suppressPackageStartupMessages(library(ArchR))
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)

GENOMEDIR='/home/bnphan/resources/genomes/rn6'

###########################################################
###### create the genome annotation files for ArchR #######
genomeAnnotation = createGenomeAnnotation(genome = BSgenome.Rnorvegicus.UCSC.rn6, 
                                          filter = TRUE, filterChr = c("chrM"))
chromSizes = genomeAnnotation$chromSizes
genome(chromSizes) <- "rn6"

# to avoid ends of chromosomes tiled regions
restrict = chromSizes
start(restrict) = start(restrict)+1e6
end(restrict) = end(restrict)-1e6

blacklist = import(file.path(GENOMEDIR,'rn6_liftOver_mm10-blacklist.v2.bed'))
blacklist = blacklist[!grepl('NW',seqnames(blacklist))]
blacklist <- sort(sortSeqlevels(blacklist), ignore.strand = TRUE)
seqlevels(blacklist) <- seqlevels(chromSizes)
seqlengths(blacklist) = seqlengths(chromSizes)
genome(blacklist) <- "rn6"
blacklist = intersect(blacklist, restrict)

genomeAnnotation$blacklist = blacklist

#################################################################
###### load in all the genomes, annotations, for rn6 #######
txdb_sqlite_fn = file.path(GENOMEDIR, 'rn6_liftoff_mm10.sqlite')
if(!file.exists(txdb_sqlite_fn)){
  # load in the gene annotation and save to sqlite 
  txdb = makeTxDbFromGFF(
    file.path(GENOMEDIR,'rn6_liftoff_mm10_RefSeq.gtf'), 
    organism = 'Rattus norvegicus', dbxrefTag = 'Dbxref') # for Rattus norvegicus
  seqlevels(txdb) <- seqlevels(chromSizes)
  saveDb(txdb, txdb_sqlite_fn)
}else {
  txdb = loadDb(txdb_sqlite_fn)
}  

txdb_rda = file.path(GENOMEDIR,'rn6_liftoff_mm10._genes_exons_TSS.rda')
if(!file.exists(txdb_rda)){
  # format genes
  library(org.Rn.eg.db)
  genes <- GenomicFeatures::genes(txdb)
  genes = genes[ !duplicated(genes)]
  genes <- dropSeqlevels(genes, grep('NW',seqlevels(genes), value = T),
                         pruning.mode="coarse")
  mcols(genes)$symbol <- mcols(genes)$gene_id
  genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
  seqlevels(genes, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(genes) = end(chromSizes)
  genome(genes) <- "rn6"
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
  genome(exons) <- "rn6"
  genes = genes[genes$symbol %in% unique(exons$symbol)]
  # exons = subsetByOverlaps(exons, restrict)
  
  # TSS
  TSS <- resize(genes, 1, "start")
  TSS <- sort(sortSeqlevels(TSS), ignore.strand = TRUE)
  TSS <- dropSeqlevels(TSS, grep('NW',seqlevels(TSS), value = T),
                       pruning.mode="coarse")
  seqlengths(TSS) = end(chromSizes)
  genome(TSS) <- "rn6"
  # TSS = subsetByOverlaps(TSS, restrict)
  
  # save the files
  save(genes, exons, TSS, file = txdb_rda)
} else {
  load(txdb_rda)
}

# make ArchR gene annotation 
geneAnnotation <- createGeneAnnotation(genes = genes, exons = exons, TSS = TSS)

save(genomeAnnotation, geneAnnotation, file = 
       file.path(GENOMEDIR,'rn6_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda'))


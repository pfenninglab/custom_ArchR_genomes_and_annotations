#!/usr/bin/env Rscript

#' Create ArchR genome and gene annotations from liftoff results
#' @param genome Target genome name (e.g. "rheMac10")
#' @param source_genome Source genome name (e.g. "hg38")
#' @param annot Annotation name (e.g. "gencode.v47.basic")
#' @param gtf_path Path to input GTF file
#' @param blacklist Optional path to blacklist BED file
#' @param outdir Output directory
#' @param tempdir Temporary directory for processing
#' @param restrict_bounds Margin size to restrict from chromosome ends (default: 1e6)

suppressPackageStartupMessages({
  library(optparse)
  library(ArchR)
  library(GenomicFeatures)
  library(rtracklayer)
  library(AnnotationDbi)
  library(tidyverse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-g", "--genome"), type="character", help="Target genome name"),
  make_option(c("-s", "--source"), type="character", help="Source genome name"),
  make_option(c("-p", "--species"), type="character", help="Species name"),
  make_option(c("-a", "--annot"), type="character", help="Annotation name"),
  make_option(c("-i", "--input"), type="character", help="Input GTF file"),
  make_option(c("-n", "--bsgenome"), type="character", help="BSgenome package name"),
  make_option(c("-b", "--blacklist"), type="character", default=NULL, help="Blacklist BED"),
  make_option(c("-m", "--chrM"), type="character", default='chrM', help="mitochrondrial chromosome name"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-t", "--temp"), type="character", help="Temp directory"),
  make_option(c("-r", "--restrict"), type="numeric", default=1e5, 
              help="Margin size to restrict from chromosome ends")
)

#' Calculate safe boundaries for gene scoring
#' @param chromSizes GRanges object with chromosome sizes
#' @param upstream Max upstream extension (default 100000)
#' @param downstream Max downstream extension (default 100000)
#' @param geneUpstream Gene body upstream extension (default 5000)
#' @param extraPadding Additional safety padding percentage (default 0.2)
#' @return List containing restricted genes and safe boundaries
#' See docs for https://www.archrproject.com/reference/addGeneScoreMatrix.html
calculateSafeBoundaries <- function(chromSizes, 
                                    upstream = 100000,
                                    downstream = 100000,
                                    geneUpstream = 5000,
                                    extraPadding = 0.2) {
  
  # Calculate maximum required padding
  maxPad <- max(c(upstream, downstream)) + geneUpstream
  safePad <- ceiling(maxPad * (1 + extraPadding))
  
  # Create restriction boundaries
  restrict <- chromSizes
  ind = width(restrict) > safePad
  start(restrict)[ind] <- start(restrict)[ind] + safePad
  end(restrict)[ind] <- end(restrict)[ind] - safePad
  
  return(list(
    restrict = restrict,
    padding = safePad
  ))
}


# Get genome name from command line args
opt <- parse_args(OptionParser(option_list=option_list))
genome_name <- opt$genome  # This will be e.g. "rheMac10", "rn7", etc.

# Construct output filenames
base_name <- paste(opt$genome, opt$source, opt$annot, sep="-")
sqlite_file <- file.path(opt$outdir, paste0(base_name, ".sqlite"))
features_file <- file.path(opt$outdir, paste0(base_name, "-genes_exons_TSS.rda"))
archr_file <- file.path(opt$outdir, paste0(base_name, "-ArchR_annotations.rda"))

# Create output directories
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(opt$temp, recursive=TRUE, showWarnings=FALSE)

# Load BSgenome package dynamically
tryCatch({
  bsgenome_pkg <- opt$bsgenome
  suppressPackageStartupMessages(require(bsgenome_pkg, character.only=TRUE))
}, error=function(e) {
  stop(paste("Failed to load BSgenome package for", opt$genome))
})

if(!is.null(opt$blacklist)) opt$blacklist = rtracklayer::import(opt$blacklist)



###########################################################
###### create the genome annotation files for ArchR #######
# Create genome annotation and get chromosome sizes
genomeAnnotation = createGenomeAnnotation(
  genome = get(bsgenome_pkg),
  filter = TRUE, 
  blacklist = opt$blacklist,
  filterChr = opt$chrM
)

chromSizes = genomeAnnotation$chromSizes
genome(chromSizes) <- opt$genome

# Process gene annotations with checkpointing
if (!file.exists(sqlite_file)) {
  message("Creating TxDb from GTF...")
  # load in the gene annotation and save to sqlite 
  txdb <- makeTxDbFromGFF(opt$input, organism=opt$species, dbxrefTag = 'Dbxref')
  seqlevels(txdb) <- seqlevels(chromSizes)
  saveDb(txdb, sqlite_file)
} else {
  message("Loading existing TxDb...")
  txdb <- loadDb(sqlite_file)
}

if (!file.exists(features_file)) {
  message("Processing gene features...")
  gtf=import.gff(opt$input, feature.type = 'transcript')
  gene_name = as.data.frame(gtf) %>% distinct(gene_id, gene_name) %>% deframe()
  
  # Process genes first
  message("Getting Genes...")
  genes <- genes(txdb)
  genes <- genes[!duplicated(genes)]
  mcols(genes)$symbol <- gene_name[mcols(genes)$gene_id] # For liftoff annotations, gene_id is the symbol
  names(genes) <- NULL
  
  # Apply chromosome filtering and genome metadata
  genes <- keepSeqlevels(genes, seqlevels(chromSizes), pruning.mode="coarse")
  seqlevels(genes) <- seqlevels(chromSizes)
  seqlengths(genes) <- end(chromSizes)
  genome(genes) <- opt$genome
  
  # Calculate safe boundaries and filter genes
  safeRegions <- calculateSafeBoundaries(
    chromSizes = chromSizes,
    upstream = 100000,
    downstream = 100000,
    geneUpstream = 5000,
    extraPadding = 0.2
  )
  
  # Apply restrictions and trim
  restrict <- safeRegions$restrict
  genes <- subsetByOverlaps(genes, restrict)
  genes <- trim(genes)
  
  # Process exons matching createGeneAnnotation style
  message("Getting Exons...")
  exons <- unlist(exonsBy(txdb, by="tx"))
  exons$exon_id <- paste0(seqnames(exons), ':', start(exons), '-', end(exons))
  
  # use the gtf to get the gene names
  gtf2=import.gff(opt$input, feature.type = 'exon')
  exon_name = gtf2 %>% as.data.frame() %>%
    mutate(exon_id = paste0(seqnames, ':', start, '-', end)) %>% 
    distinct(exon_id, gene_name) %>% deframe()
  
  mcols(exons)$symbol <- exon_name[mcols(exons)$exon_id]
  
  # Clean up exon metadata
  names(exons) <- NULL
  mcols(exons)$exon_id <- NULL
  mcols(exons)$exon_name <- NULL
  mcols(exons)$exon_rank <- NULL
  mcols(exons)$tx_id <- NULL
  
  # Filter and sort exons
  exons <- sort(sortSeqlevels(exons))
  exons <- keepSeqlevels(exons, seqlevels(chromSizes), pruning.mode="coarse")
  seqlevels(exons) <- seqlevels(chromSizes)
  seqlengths(exons) <- end(chromSizes)
  genome(exons) <- opt$genome
  
  # Only keep exons with valid genes
  exons <- exons[!is.na(exons$symbol) & 
                   !duplicated(exons) & 
                   exons$symbol %in% genes$symbol]
  
  # Apply restrictions and trim
  exons <- subsetByOverlaps(exons, restrict)
  exons <- trim(exons)
  
  # Update genes based on exon symbols
  genes <- genes[genes$symbol %in% unique(exons$symbol)]
  
  # Create TSS from filtered genes
  message("Getting TSS...")
  TSS <- resize(genes, width=1, fix="start")  
  TSS <- sort(sortSeqlevels(TSS))
  TSS <- trim(TSS)
  
  # Save processed features
  save(genes, exons, TSS, restrict, file=features_file)
} else {
  message("Loading existing gene features...")
  load(features_file)
}

# Process blacklist if provided
if (!is.null(opt$blacklist)) {
  message("Processing blacklist...")
  blacklist = genomeAnnotation$blacklist
  seqlevels(blacklist) <- seqlevels(chromSizes)
  seqlengths(blacklist) <- seqlengths(chromSizes)
  blacklist <- sort(sortSeqlevels(blacklist))
  genome(blacklist) <- opt$genome
  blacklist <- trim(blacklist)
  blacklist <- GenomicRanges::intersect(blacklist, restrict)
  genomeAnnotation$blacklist <- blacklist
}

# Create ArchR gene annotation and save final output
message("Creating ArchR annotations...")
geneAnnotation <- createGeneAnnotation(genes=genes, exons=exons, TSS=TSS)
save(genomeAnnotation, geneAnnotation, file=archr_file)

message("Done!")
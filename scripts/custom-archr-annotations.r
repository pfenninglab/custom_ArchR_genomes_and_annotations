#!/usr/bin/env Rscript

# Script to generate ArchR genome and gene annotations from lifted annotations
# Author: Claude
# Usage: Rscript generate_archr_annotations.R --genome [genome] --gtf [gtf_path] --blacklist [blacklist_path] --outdir [output_dir]

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ArchR)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(futile.logger)
})

# Setup logging
setup_logging <- function(outdir) {
  log_file <- file.path(outdir, "archr_annotations.log")
  flog.appender(appender.file(log_file))
  flog.info("Starting ArchR annotation generation")
}

# Parse command line arguments
option_list <- list(
  make_option("--genome", type="character", help="Genome name (e.g. mCalJac1)"),
  make_option("--gtf", type="character", help="Path to GTF annotation file"),
  make_option("--blacklist", type="character", default=NULL, help="Path to blacklist regions"),
  make_option("--outdir", type="character", help="Output directory"),
  make_option("--bsgenome", type="character", help="BSgenome package name")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate required arguments
required <- c("genome", "gtf", "outdir", "bsgenome")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse=", "))
}

# Create output directory if it doesn't exist
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
setup_logging(opt$outdir)

# Function to load BSgenome package
load_bsgenome <- function(package_name) {
  tryCatch({
    suppressPackageStartupMessages(require(package_name, character.only=TRUE))
    flog.info("Loaded BSgenome package: %s", package_name)
    get(package_name)
  }, error=function(e) {
    flog.error("Failed to load BSgenome package: %s", e$message)
    stop("BSgenome package loading failed")
  })
}

# Function to create genome annotation
create_genome_annot <- function(bsgenome, blacklist_path=NULL) {
  flog.info("Creating genome annotation")
  
  # Create basic genome annotation
  genome_annot <- tryCatch({
    createGenomeAnnotation(
      genome = bsgenome,
      filter = TRUE
    )
  }, error=function(e) {
    flog.error("Failed to create genome annotation: %s", e$message)
    stop("Genome annotation creation failed")
  })
  
  # Add blacklist if provided
  if (!is.null(blacklist_path) && file.exists(blacklist_path)) {
    flog.info("Adding blacklist regions from: %s", blacklist_path)
    blacklist <- tryCatch({
      bl <- import(blacklist_path)
      bl <- bl[!grepl('NW|NT|GL|KI|JH', seqnames(bl))]
      seqlevels(bl) <- seqlevels(genome_annot$chromSizes)
      seqlengths(bl) <- end(genome_annot$chromSizes)
      genome(bl) <- opt$genome
      bl
    }, error=function(e) {
      flog.error("Failed to process blacklist: %s", e$message)
      NULL
    })
    
    if (!is.null(blacklist)) {
      genome_annot$blacklist <- blacklist
    }
  }
  
  return(genome_annot)
}

# Function to create gene annotation
create_gene_annot <- function(gtf_path, genome_annot) {
  flog.info("Creating gene annotation from GTF: %s", gtf_path)
  
  # Create TxDb from GTF
  txdb <- tryCatch({
    makeTxDbFromGFF(gtf_path, organism=opt$genome)
  }, error=function(e) {
    flog.error("Failed to create TxDb from GTF: %s", e$message)
    stop("TxDb creation failed")
  })
  
  # Extract genes and filter those near chromosome ends
  genes <- tryCatch({
    genes <- genes(txdb)
    genes <- genes[!duplicated(genes)]
    
    # Set seqlevels and lengths
    seqlevels(genes) <- seqlevels(genome_annot$chromSizes)
    seqlengths(genes) <- end(genome_annot$chromSizes)
    genome(genes) <- opt$genome
    
    # Filter genes near chromosome ends (120kb from either end)
    BUFFER_SIZE <- 120000
    chrom_sizes <- end(genome_annot$chromSizes)
    
    # Create GRanges objects for chromosome start and end regions to exclude
    chrom_starts <- GRanges(
      seqnames = names(chrom_sizes),
      ranges = IRanges(start = 1, end = BUFFER_SIZE)
    )
    chrom_ends <- GRanges(
      seqnames = names(chrom_sizes),
      ranges = IRanges(
        start = chrom_sizes - BUFFER_SIZE,
        end = chrom_sizes
      )
    )
    
    # Find genes that overlap with excluded regions
    genes_near_start <- overlapsAny(genes, chrom_starts)
    genes_near_end <- overlapsAny(genes, chrom_ends)
    
    # Keep only genes that are not near chromosome ends
    filtered_genes <- genes[!(genes_near_start | genes_near_end)]
    
    flog.info("Filtered out %d genes within %d bp of chromosome ends",
              length(genes) - length(filtered_genes),
              BUFFER_SIZE)
    
    filtered_genes
  }, error=function(e) {
    flog.error("Failed to process genes: %s", e$message)
    stop("Gene processing failed")
  })
  
  # Extract exons and filter based on filtered genes
  exons <- tryCatch({
    exons <- unlist(exonsBy(txdb, by="tx"))
    exons$symbol <- names(exons)
    exons <- exons[!duplicated(exons)]
    seqlevels(exons) <- seqlevels(genome_annot$chromSizes)
    seqlengths(exons) <- end(genome_annot$chromSizes)
    genome(exons) <- opt$genome
    
    # Only keep exons that belong to our filtered genes
    valid_gene_ids <- mcols(filtered_genes)$gene_id
    exons <- exons[exons$symbol %in% valid_gene_ids]
    
    flog.info("Kept %d exons associated with filtered genes", length(exons))
    exons
  }, error=function(e) {
    flog.error("Failed to process exons: %s", e$message)
    stop("Exon processing failed") 
  })
  
  # Create TSS from filtered genes
  tss <- tryCatch({
    # Create TSS only from our filtered genes
    tss <- resize(filtered_genes, 1, "start")
    seqlevels(tss) <- seqlevels(genome_annot$chromSizes)
    seqlengths(tss) <- end(genome_annot$chromSizes)
    genome(tss) <- opt$genome
    
    flog.info("Created %d TSS sites from filtered genes", length(tss))
    tss
  }, error=function(e) {
    flog.error("Failed to create TSS: %s", e$message)
    stop("TSS creation failed")
  })
  
  # Create gene annotation
  gene_annot <- createGeneAnnotation(genes=genes, exons=exons, TSS=tss)
  
  return(gene_annot)
}

# Main execution
main <- function() {
  # Load BSgenome
  bsgenome <- load_bsgenome(opt$bsgenome)
  
  # Create genome annotation
  genome_annot <- create_genome_annot(bsgenome, opt$blacklist)
  
  # Create gene annotation
  gene_annot <- create_gene_annot(opt$gtf, genome_annot)
  
  # Save annotations
  output_file <- file.path(opt$outdir, paste0(opt$genome, "_ArchR_annotations.rda"))
  flog.info("Saving annotations to: %s", output_file)
  save(genome_annot, gene_annot, file=output_file)
  
  flog.info("Successfully completed ArchR annotation generation")
}

# Execute main function with error handling
tryCatch({
  main()
}, error=function(e) {
  flog.error("Fatal error: %s", e$message)
  quit(status=1)
})
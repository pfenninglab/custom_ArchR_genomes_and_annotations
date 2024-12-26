#!/usr/bin/env Rscript
#
# Simple script to create BSgenome package from 2bit file
# 
# Usage:
#   Rscript script.R <proj_dir> <species_name> <genome_id> <twobit_path> <output_dir> <pkg_name>
#
# Example:
#   Rscript script.R /path/to/project "Macaca nemestrina" mMacNem1 /path/to/genome.2bit /output/dir BSgenome.pkg.name
#
# Author: BaDoi Phan
# Date: December 2024

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if (length(args) != 6) {
  stop("Usage: Rscript script.R <proj_dir> <species_name> <genome_id> <twobit_path> <output_dir> <pkg_name>")
}

# Store arguments in named list for clarity
params <- list(
  proj_dir = args[1],    # Project root directory
  species_name = args[2], # Scientific name (e.g. "Macaca nemestrina")
  genome_id = args[3],    # Genome assembly ID (e.g. "mMacNem1") 
  twobit_path = args[4],  # Path to .2bit file
  output_dir = args[5],   # Output directory for package
  pkg_name = args[6]      # Name for BSgenome package
)

# Source the patched functions
patch_file <- file.path(params$proj_dir, 'scripts', 'bsgenomeforge-patch.R')
if (!file.exists(patch_file)) {
  stop("Patch file not found at: ", patch_file)
}
source(patch_file)

# Print debug info
cat("\nDEBUG INFO:\n")
cat("Project directory:", params$proj_dir, "\n")
cat("Species name:", params$species_name, "\n")
cat("Genome ID:", params$genome_id, "\n")
cat("2bit file exists:", file.exists(params$twobit_path), "\n")
cat("Output dir exists:", dir.exists(params$output_dir), "\n")

# Create output directory if it doesn't exist
if (!dir.exists(params$output_dir)) {
  message("Creating output directory:", params$output_dir)
  dir.create(params$output_dir, recursive = TRUE)
  
  dir.create(file.path(params$output_dir, params$pkg_name, "inst", "extdata"), 
             recursive = TRUE, showWarnings = FALSE)
}

# Create BSgenome package
message("\nCreating BSgenome package...")
forgeBSgenomeDataPkgFromTwobitFile(
  # Required file path parameters
  filepath = params$twobit_path,
  organism = params$species_name,
  provider = "Custom",
  genome = params$genome_id,
  
  # Package metadata
  pkg_maintainer = "BaDoi Phan <badoi@me.com>",
  pkg_license = "Artistic-2.0",
  pkg_author = 'BaDoi Phan',
  circ_seqs = character(0),
  pkg_version = "1.0.0",
  
  # Output location
  destdir = params$output_dir
)

message("Done!")
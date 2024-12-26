# Check and install required packages with specific versions
if (!require("BSgenomeForge", quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  message("Installing BSgenomeForge...")
  # Install specific release versions to ensure compatibility 
  remotes::install_github("Bioconductor/BSgenomeForge@RELEASE_3_20")
  remotes::install_github("Bioconductor/Biobase@RELEASE_3_20")
}

# Load required packages silently
suppressPackageStartupMessages({
  library(BSgenomeForge)
  library(Biobase)
})

# Main forge function - creates BSgenome data package from 2bit file
forgeBSgenomeDataPkgFromTwobitFile <- function(filepath,
                                               organism, provider, genome,
                                               pkg_maintainer, pkg_author,
                                               pkg_version="1.0.0",
                                               pkg_license="Artistic-2.0",
                                               seqnames=NULL,
                                               circ_seqs=character(0),
                                               destdir=".") {
  # Input validation
  if (!isSingleString(organism) || organism == "")
    stop(wmsg("'organism' must be a single (non-empty) string"))
  # ... other validation checks ...
  
  # Get sequence names from 2bit file and validate
  file_seqnames <- BSgenomeForge:::.get_seqnames_from_twobit_file(filepath)
  BSgenomeForge:::.check_seqnames(seqnames, filepath, file_seqnames)
  circ_seqs <- BSgenomeForge:::.normarg_circ_seqs(circ_seqs, filepath, file_seqnames,
                                                  seqnames=seqnames)
  
  # Sort sequences if needed
  if (!is.null(seqnames)) {
    sorted_twobit_file <- file.path(tempdir(), "single_sequences.2bit")
    filepath <- BSgenomeForge:::.sort_twobit_file1(filepath, sorted_twobit_file, seqnames)
  }
  
  # Format package metadata
  organism <- BSgenomeForge:::format_organism(organism)
  abbr_organism <- BSgenomeForge:::abbreviate_organism_name(organism)
  pkgname <- BSgenomeForge:::make_pkgname(abbr_organism, provider, genome)
  pkg_title <- BSgenomeForge:::.make_pkgtitle(organism, provider, genome)
  pkg_desc <- BSgenomeForge:::.make_pkgdesc(organism, provider, genome)
  biocview <- BSgenomeForge:::organism2biocview(organism)
  forge_function <- "forgeBSgenomeDataPkgFromTwobitFile"
  
  # Create the package
  create_2bit_BSgenome_datapkg(filepath, pkgname,
                               BSgenome_objname=abbr_organism,
                               pkg_title=pkg_title,
                               pkg_desc=pkg_desc,
                               pkg_version=pkg_version,
                               pkg_author=pkg_author,
                               pkg_maintainer=pkg_maintainer,
                               pkg_license=pkg_license,
                               organism=organism,
                               provider=provider,
                               genome=genome,
                               organism_biocview=biocview,
                               forge_function=forge_function,
                               circ_seqs=circ_seqs,
                               destdir=destdir,
                               move_twobit_file=!is.null(seqnames))
}

# Helper function to create the actual package structure
create_2bit_BSgenome_datapkg <- function(twobit_path, pkgname, BSgenome_objname,
                                         pkg_title, pkg_desc, pkg_version,
                                         pkg_author, pkg_maintainer, pkg_license,
                                         organism, provider, genome, organism_biocview, 
                                         forge_function, circ_seqs, destdir=".", 
                                         move_twobit_file=FALSE) {
  # Validate inputs
  stopifnot(isSingleString(twobit_path),
            isSingleString(pkgname),
            isSingleString(destdir),
            isTRUEorFALSE(move_twobit_file))
  
  # Get package template directory
  origdir <- system.file("pkgtemplates", "2bit_BSgenome_datapkg",
                         package="BSgenomeForge", mustWork=TRUE)
  
  # Create package values list
  symValues <- list(BSGENOMEOBJNAME=BSgenome_objname,
                    PKGTITLE=pkg_title,
                    PKGDESCRIPTION=pkg_desc,
                    PKGVERSION=pkg_version,
                    PKGAUTHOR=pkg_author,
                    PKGMAINTAINER=pkg_maintainer,
                    PKGLICENSE=pkg_license,
                    ORGANISM=organism,
                    PROVIDER=provider,
                    GENOME=genome,
                    ORGANISMBIOCVIEW=organism_biocview,
                    CIRCSEQS=BSgenomeForge:::build_Rexpr_as_string(circ_seqs),
                    FORGEFUN=forge_function)
  
  # Validate all values are strings
  stopifnot(all(vapply(symValues, isSingleString, logical(1))))
  
  # Create the package structure
  pkg_dir <- Biobase:::createPackage(pkgname, destdir, origdir, symValues,
                                     unlink=TRUE, quiet=FALSE)[[1L]]
  
  # Fix: Create extdata directory structure explicitly
  dir.create(file.path(pkg_dir, "inst", "extdata"), recursive=T, showWarnings=F)
  
  # Copy or move the 2bit file
  to <- file.path(pkg_dir, "inst", "extdata", "single_sequences.2bit")
  if (move_twobit_file) {
    file.rename(twobit_path, to)
  } else {
    stopifnot(file.copy(twobit_path, to))
  }
  
  invisible(pkg_dir)
}
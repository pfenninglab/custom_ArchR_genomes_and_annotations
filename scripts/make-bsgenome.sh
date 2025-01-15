#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL
#SBATCH --mem=24G
#SBATCH --error=logs/%x_%j.txt
#SBATCH --output=logs/%x_%j.txt

set -euo pipefail

# Function to print usage
usage() {
    echo "Usage: $0 -s SPECIES -g GENOME -f FASTA -o OUTPUT_DIR -n PKG_NAME -p PROJECT_DIR [-w SCRATCH_DIR]" 1>&2
    echo ""
    echo "Required arguments:"
    echo "  -s SPECIES      Species name (e.g., 'Macaca nemestrina')"
    echo "  -g GENOME       Genome version (e.g., 'mMacNem1')"
    echo "  -f FASTA       Path to genome FASTA file"
    echo "  -o OUTPUT_DIR   Output directory for BSgenome package"
    echo "  -n PKG_NAME     Package name (must start with 'BSgenome.')"
    echo "  -p PROJECT_DIR  Project directory containing scripts/"
    echo ""
    echo "Optional arguments:"
    echo "  -w SCRATCH_DIR  Scratch directory (default: /scratch/\$USER)"
    echo "  -h             Show this help message"
    exit 1
}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to validate species name format
validate_species_name() {
    local species="$1"
    if [[ ! "$species" =~ ^[A-Z][a-z]+[[:space:]][a-z]+$ ]]; then
        log "ERROR: Species name must be in format 'Genus species' (e.g., 'Macaca nemestrina')"
        exit 1
    fi
}

# Function to validate package name format
validate_package_name() {
    local pkg_name="$1"
    if [[ ! "$pkg_name" =~ ^BSgenome\. ]]; then
        log "ERROR: Package name must start with 'BSgenome.'"
        exit 1
    fi
}

make_names() {
    local input="$1"
    
    # Step 1: Replace non-alphanumeric characters (except underscores) with a period
    local name=$(echo "$input" | sed 's/[^a-zA-Z0-9.]//g')

    # Step 2: Replace multiple consecutive periods with a single period
    name=$(echo "$name" | sed 's/\.\+/\./g')

    # Step 3: Trim leading and trailing periods
    name=$(echo "$name" | sed 's/^\.//; s/\.$//')

    # Step 4: Ensure the name is not empty; use "X" as a fallback
    if [[ -z "$name" ]]; then
        name="X"
    fi

    echo "$name"
}

# Parse command line arguments
while getopts "s:g:f:o:n:p:w:h" opt; do
    case $opt in
        s) SPECIES="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        f) FASTA="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        n) PKG_NAME="$OPTARG" ;;
        p) PROJECT_DIR="$OPTARG" ;;
        w) SCRATCH_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Activate conda environment with liftoff
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate custom_genes

# Log all variables
log "SPECIES: $SPECIES"
log "GENOME: $GENOME"
log "FASTA: $FASTA"
log "OUTPUT_DIR: $OUTPUT_DIR"
log "PKG_NAME: $PKG_NAME"
log "PROJECT_DIR: $PROJECT_DIR"
log "SCRATCH_DIR: ${SCRATCH_DIR:-/scratch/$USER}"

# Check required arguments
if [[ -z "${SPECIES:-}" || -z "${GENOME:-}" || -z "${FASTA:-}" || \
      -z "${OUTPUT_DIR:-}" || -z "${PKG_NAME:-}" || -z "${PROJECT_DIR:-}" ]]; then
    usage
fi

# Validate inputs
validate_species_name "$SPECIES"
validate_package_name "$PKG_NAME"

# Check if FASTA file exists and is readable
if [[ ! -f "$FASTA" ]]; then
    log "ERROR: FASTA file not found at $FASTA"
    exit 1
fi

if [[ ! -r "$FASTA" ]]; then
    log "ERROR: FASTA file exists but is not readable at $FASTA"
    exit 1
fi

# Check if package is already installed
if Rscript -e "if('${PKG_NAME}' %in% rownames(installed.packages())) quit(status = 0) else quit(status = 1)"; then
    log "Package ${PKG_NAME} is already installed, skipping all steps..."
    exit 0
fi

# Check if tar.gz already exists
TAR_FILE="${OUTPUT_DIR}/${PKG_NAME}_1.0.0.tar.gz"
if [[ -f "$TAR_FILE" ]]; then
    log "Package tar.gz already exists at ${TAR_FILE}, installing..."
    R CMD INSTALL "$TAR_FILE"
    if [[ $? -eq 0 ]]; then
        log "Package installed successfully"
        exit 0
    else
        log "Package installation failed, continuing with full package creation..."
    fi
fi

# Set default scratch directory if not specified
SCRATCH_DIR=${SCRATCH_DIR:-/scratch/$USER}
SCRIPTS_DIR="$PROJECT_DIR/scripts"

# Verify project directory structure
if [[ ! -d "$SCRIPTS_DIR" ]]; then
    log "ERROR: Scripts directory not found at $SCRIPTS_DIR"
    exit 1
fi

# Set up working directories
TEMP_DIR="${SCRATCH_DIR}/${GENOME}-bsgenome-${SLURM_JOB_ID}"
mkdir -p "$TEMP_DIR" "$OUTPUT_DIR"

cd $TEMP_DIR

# Convert FASTA to 2bit
TWOBIT_FILE="$PROJECT_DIR/output/genomes/$GENOME/$GENOME.2bit"
if [[ ! -f "$TWOBIT_FILE" ]]; then
    log "Converting FASTA to 2bit format..."

    # Copy FASTA to temporary directory
    if ! rsync "$FASTA" "$TEMP_DIR/"; then
        log "ERROR: Failed to copy FASTA file to temporary directory"
        rm -rf "$TEMP_DIR"
        exit 1
    fi
    FASTA=$(basename "$FASTA")

    if [[ "$FASTA" =~ \.gz$ ]]; then
        gunzip "$TEMP_DIR/$FASTA"
        FASTA=$(basename "$FASTA" .gz)
    fi

    # Convert
    TWOBIT_FILE=$(basename "$FASTA")
    faToTwoBit "$FASTA" "$TWOBIT_FILE"
    rsync -Paq $TWOBIT_FILE $PROJECT_DIR/output/genomes/${GENOME}/
else 
    log "2bit file already exists, grabbing from output directory..."
    rsync -Paq $PROJECT_DIR/output/genomes/$GENOME/$GENOME.2bit $TEMP_DIR
    TWOBIT_FILE="$TEMP_DIR/$GENOME.2bit"
fi

# Create BSgenome package using R script
log "Creating BSgenome package..."
Rscript "$SCRIPTS_DIR/make-bsgenome.R" \
    "$PROJECT_DIR" \
    "$SPECIES" \
    "$GENOME" \
    "$TWOBIT_FILE" \
    "$TEMP_DIR" \
    "$PKG_NAME"

PKG_NAME=$(make_names $PKG_NAME)

# Build and install R package
if [[ -d "$PKG_NAME" ]]; then
    log "Building R package..."    
    R CMD build "$PKG_NAME"
    R CMD INSTALL "${PKG_NAME}_1.0.0.tar.gz"

    # Copy package tar to output directory
    rsync -Pav "${PKG_NAME}_1.0.0.tar.gz" ${OUTPUT_DIR}/
fi

# Cleanup
log "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

log "BSgenome package creation complete for $GENOME"
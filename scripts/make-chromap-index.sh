#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --time=4:00:00
#SBATCH --export=ALL
#SBATCH --mem=24G
#SBATCH --error=logs/%x_%j.txt
#SBATCH --output=logs/%x_%j.txt

set -euo pipefail

# Function to print usage
usage() {
    echo "Usage: $0 -g GENOME -f FASTA [-w SCRATCH_DIR] [-o OUTPUT_DIR]" 1>&2
    echo ""
    echo "Build Chromap index for snATAC-seq analysis"
    echo ""
    echo "Required arguments:"
    echo "  -g GENOME       Genome version (e.g., 'mm10', 'hg38', 'rn6')"
    echo "  -f FASTA       Path to genome FASTA file"
    echo ""
    echo "Optional arguments:"
    echo "  -o OUTPUT_DIR   Output directory (default: /path/to/genomes/chromap/GENOME)"
    echo "  -w SCRATCH_DIR  Scratch directory (default: /scratch/\$USER)"
    echo "  -h             Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -g mm10 -f /path/to/genomes/mm10/mm10.fa"
    exit 1
}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Parse command line arguments
while getopts "g:f:o:w:h" opt; do
    case $opt in
        g) GENOME="$OPTARG" ;;
        f) FASTA="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        w) SCRATCH_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [[ -z "${GENOME:-}" || -z "${FASTA:-}" ]]; then
    usage
fi

# Set default directories if not specified
SCRATCH_DIR=${SCRATCH_DIR:-/scratch/$USER}
OUTPUT_DIR=${OUTPUT_DIR:-/home/$USER/resources/genomes/$GENOME/chromap}

# Log all variables
log "GENOME: $GENOME"
log "FASTA: $FASTA"
log "OUTPUT_DIR: $OUTPUT_DIR"
log "SCRATCH_DIR: $SCRATCH_DIR"

# Check if FASTA file exists and is readable
if [[ ! -f "$FASTA" ]]; then
    log "ERROR: FASTA file not found at $FASTA"
    exit 1
fi

if [[ ! -r "$FASTA" ]]; then
    log "ERROR: FASTA file exists but is not readable at $FASTA"
    exit 1
fi

# Check if output directory exists and create if needed
mkdir -p "$OUTPUT_DIR"

# Check if index already exists
if [[ -f "${OUTPUT_DIR}/index" ]]; then
    log "Chromap index already exists at ${OUTPUT_DIR}/index"
    exit 0
fi

# Set up working directories
TEMP_DIR="${SCRATCH_DIR}/${GENOME}-chromap-${SLURM_JOB_ID}"
mkdir -p "$TEMP_DIR"

# Activate required environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate custom_genes

# Copy FASTA to temporary directory
log "Copying FASTA to temporary directory..."
if ! rsync -P "$FASTA" "$TEMP_DIR/"; then
    log "ERROR: Failed to copy FASTA file to temporary directory"
    rm -rf "$TEMP_DIR"
    exit 1
fi

cd "$TEMP_DIR"

# Handle compressed FASTA
FASTA_BASE=$(basename "$FASTA")
if [[ "$FASTA_BASE" =~ \.gz$ ]]; then
    log "Decompressing FASTA file..."
    gunzip "$FASTA_BASE"
    FASTA_BASE="${FASTA_BASE%.gz}"
fi

# Create Chromap index for snATAC-seq
log "Creating Chromap index for snATAC-seq..."
chromap -i -r "$FASTA_BASE" -o index

# Copy results back to output directory
log "Copying index to output directory..."
rsync -Paq --remove-source-files index $OUTPUT_DIR

# Cleanup
log "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

log "Chromap index creation complete for $GENOME"

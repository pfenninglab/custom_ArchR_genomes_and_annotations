#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --time=4:00:00
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --error=logs/%x_%j.txt
#SBATCH --output=logs/%x_%j.txt

set -euo pipefail

# Function to print usage
usage() {
    echo "Usage: $0 -g GENOME -s SOURCE_GENOME -a ANNOTATION [-p PROJECT_DIR] [-w SCRATCH_DIR] [-c CONFIG]" 1>&2
    echo ""
    echo "Create ArchR genome and gene annotation files from liftoff results"
    echo ""
    echo "Required arguments:"
    echo "  -g GENOME          Target genome name (e.g. rheMac10)"
    echo "  -s SOURCE_GENOME   Source genome name (e.g. hg38)" 
    echo "  -a ANNOTATION      Annotation name (e.g. gencode.v47.basic)"
    echo ""
    echo "Optional arguments:"
    echo "  -p PROJECT_DIR     Project directory (default: current directory)"
    echo "  -w SCRATCH_DIR     Scratch directory (default: /scratch/\$USER)"
    echo "  -c CONFIG          Path to config file (default: \$PROJECT_DIR/inst/config/target_genomes.tsv)"
    echo "  -h                 Show this help message"
    exit 1
}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to clean up on exit
cleanup() {
    if [[ -d "$TEMP_DIR" ]]; then
        log "Cleaning up temporary directory"
        rm -rf "$TEMP_DIR"
    fi
}

# Function to get genome info from config
get_genome_info() {
    local genome="$1"
    local config_file="$2"
    
    if [[ ! -f "$config_file" ]]; then
        log "ERROR: Config file not found: $config_file"
        return 1
    fi
    
    # Skip header line and find matching genome
    # Expected columns: common_name genome_name species_name ... bsgenome chrM
    awk -F'\t' -v genome="$genome" '
        NR == 1 {
            for (i=1; i<=NF; i++) {
                if ($i == "species_name") species_col = i
                if ($i == "bsgenome") bsgenome_col = i
                if ($i == "chrM") chrm_col = i
            }
            if (!species_col || !bsgenome_col) {
                print "ERROR: Required columns missing from config" > "/dev/stderr"
                exit 1
            }
        }
        NR > 1 && $2 == genome {
            print $species_col "\t" $bsgenome_col "\t" (chrm_col ? $chrm_col : "chrM")
        }' "$config_file"
}


# Set cleanup trap
trap cleanup EXIT

# Parse command line arguments
while getopts "g:s:a:p:w:c:h" opt; do
    case $opt in
        g) GENOME="$OPTARG" ;;
        s) SOURCE_GENOME="$OPTARG" ;;
        a) ANNOTATION="$OPTARG" ;;
        p) PROJECT_DIR="$OPTARG" ;;
        w) SCRATCH_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [[ -z "${GENOME:-}" || -z "${SOURCE_GENOME:-}" || -z "${ANNOTATION:-}" ]]; then
    usage
fi

# Set default directories if not specified
PROJECT_DIR=${PROJECT_DIR:-$(pwd)}
SCRATCH_DIR=${SCRATCH_DIR:-/scratch/$USER}
CONFIG_FILE=${CONFIG_FILE:-"$PROJECT_DIR/config/target_genomes.tsv"}

log "Using config file: $CONFIG_FILE"

# Get genome info from config
GENOME_INFO=$(get_genome_info "$GENOME" "$CONFIG_FILE")

if [[ -z "$GENOME_INFO" ]]; then
    log "ERROR: Could not find genome info for $GENOME in config file"
    exit 1
fi

# Parse genome info
IFS=$'\t' read -r species_name bsgenome CHRM <<< "$GENOME_INFO"
CHRM=${CHRM:-chrM} # Default to chrM if not specified

log "Found genome metadata:"
log "  Species: $SPECIES"
log "  BSgenome package: $BSGENOME"
log "  Mitochondrial chromosome: $CHRM"

# Set up directory structure
TEMP_DIR="${SCRATCH_DIR}/${GENOME}-archr-${SLURM_JOB_ID}"
GENOME_DIR="$PROJECT_DIR/genomes/$GENOME"

log "Creating directories..."
mkdir -p "$TEMP_DIR"/{input,temp,output} \
        "$GENOME_DIR/annotations/sqlite" \
        "$GENOME_DIR/archr"

# Construct file paths
BASE_NAME="${GENOME}-${SOURCE_GENOME}-${ANNOTATION}"
GTF_PATH="$GENOME_DIR/annotations/${BASE_NAME}.gtf.gz"
BLACKLIST_PATH="$GENOME_DIR/${GENOME}_blacklist.bed.gz"

cd "$TEMP_DIR"

# Copy and decompress input files
log "Preparing input files..."
if [[ -f "$GTF_PATH" ]]; then
    cp "$GTF_PATH" "$TEMP_DIR/input/"
    gunzip -f "$TEMP_DIR/input/$(basename "$GTF_PATH")"
else
    log "ERROR: GTF file not found at $GTF_PATH"
    exit 1
fi

# Check for blacklist and prepare if exists
BLACKLIST_ARG=""
if [[ -f "$BLACKLIST_PATH" ]]; then
    log "Found blacklist, including in processing..."
    cp "$BLACKLIST_PATH" "$TEMP_DIR/input/"
    gunzip -f "$TEMP_DIR/input/$(basename "$BLACKLIST_PATH")"
    BLACKLIST_ARG="--blacklist $TEMP_DIR/input/$(basename "${BLACKLIST_PATH%.gz}")"
fi

# Execute R script
log "Running R script to create ArchR annotations..."
Rscript "$PROJECT_DIR/R/create-archr-genome.R" \
    --genome "$GENOME" \
    --source "$SOURCE_GENOME" \
    --species "$SPECIES" \
    --annot "$ANNOTATION" \
    --input "$TEMP_DIR/input/${BASE_NAME}.gtf" \
    --bsgenome "$BSGENOME" \
    --chrM "$CHRM" \
    --outdir "$TEMP_DIR" \
    --temp "$TEMP_DIR/temp" \
    $BLACKLIST_ARG

# Copy results back to project directory
log "Copying results to project directory..."
for ext in "sqlite" "genes_exons_TSS.rda" "ArchR_annotations.rda"; do
    src="$TEMP_DIR/${BASE_NAME}-${ext}"
    if [[ -f "$src" ]]; then
        if [[ "$ext" == "sqlite" ]]; then
            dst="$GENOME_DIR/annotations/sqlite/${BASE_NAME}.${ext}"
        elif [[ "$ext" == *"ArchR"* ]]; then
            dst="$GENOME_DIR/archr/${BASE_NAME}-${ext}"
        else
            dst="$GENOME_DIR/annotations/${BASE_NAME}-${ext}"
        fi
        rsync -Paq "$src" "$dst"
    else
        log "WARNING: Expected output file not found: $(basename "$src")"
    fi
done

log "ArchR genome creation complete for $GENOME"
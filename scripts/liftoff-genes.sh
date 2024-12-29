#!/bin/bash
#SBATCH --partition=pool3-bigmem,pfen3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --error=logs/%x%j.txt
#SBATCH --output=logs/%x%j.txt

set -euo pipefail

# Function to print usage
usage() {
    echo "Usage: $0 -s SOURCE_GENOME -t TARGET_GENOME -a SOURCE_ANNOT [-p PROJECT_DIR] [-w SCRATCH_DIR]" 1>&2
    echo "Maps gene annotations between genomes using liftoff"
    echo ""
    echo "Required arguments:"
    echo "  -s SOURCE_GENOME   Source genome name (e.g. hg38)" 
    echo "  -t TARGET_GENOME   Target genome name (e.g. mm10)"
    echo "  -a SOURCE_ANNOT    Source annotation name (e.g. gencode.v39)"
    echo ""
    echo "Optional arguments:"
    echo "  -p PROJECT_DIR     Project directory (default: current directory)"
    echo "  -w SCRATCH_DIR     Scratch directory (default: /scratch/$USER)"
    exit 1
}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check for gene features in GTF
check_gene_feature() {
    local gtf_file="$1"
    log "Checking first 1000 lines for gene features in GTF..."
    
    # Check if file is compressed and count gene features
    if [[ "$gtf_file" =~ \.gz$ ]]; then
        gene_count=$(zcat "$gtf_file" | head -n 1000 | awk '$3=="gene"' | wc -l)
    else
        gene_count=$(head -n 1000 "$gtf_file" | awk '$3=="gene"' | wc -l)
    fi
        
    # Log the finding
    if [ "$gene_count" -gt 0 ]; then
        log "Found $gene_count gene features in GTF"
    else
        log "No gene features found in first 1000 lines, will use -infer_genes"
    fi

    echo "$gene_count"  # Return the actual count
}

# Function to handle input files
handle_input_file() {
    local file="$1"
    local base=$(basename "$file" .gz)
    
    if [[ -f "${file}.gz" ]]; then
        log "Using compressed ${base}.gz"
        rsync -Paq "${file}.gz" "./${base}.gz"
        gunzip "./${base}.gz"
    else
        log "Using uncompressed ${base}"
        rsync -Paq "$file" "./${base}"
    fi
}

# Parse command line arguments
while getopts "s:t:a:p:w:h" opt; do
    case $opt in
        s) SOURCE_GENOME="$OPTARG" ;;
        t) TARGET_GENOME="$OPTARG" ;;
        a) SOURCE_ANNOT="$OPTARG" ;;
        p) PROJECT_DIR="$OPTARG" ;;
        w) SCRATCH_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [[ -z "${SOURCE_GENOME:-}" ]] || [[ -z "${TARGET_GENOME:-}" ]] || [[ -z "${SOURCE_ANNOT:-}" ]]; then
    usage
fi

# requested mem in MB convert to GB and estimate cores
(( MEM_GB = ${SLURM_MEM_PER_NODE} / 1024 ))
(( N_CORES = ${MEM_GB} / 8 ))
echo -e "Memory requested for job is ${MEM_GB}GB."
echo -e "We're going to set parallelism to ${N_CORES} cores."

# Set default directories if not specified
PROJECT_DIR=${PROJECT_DIR:-$(pwd)}
SCRATCH_DIR=${SCRATCH_DIR:-/scratch/$USER}

# Set up directory structure
TEMP_DIR="$SCRATCH_DIR/${SOURCE_GENOME}-${TARGET_GENOME}-${SOURCE_ANNOT}-${SLURM_JOB_ID}"
mkdir -p "$TEMP_DIR"

# Define input/output paths
SOURCE_FA="$PROJECT_DIR/output/genomes/$SOURCE_GENOME/$SOURCE_GENOME.fa"
TARGET_FA="$PROJECT_DIR/output/genomes/$TARGET_GENOME/$TARGET_GENOME.fa"
SOURCE_GTF="$PROJECT_DIR/output/genomes/$SOURCE_GENOME/annotations/${SOURCE_ANNOT}.gtf"
OUT_DIR="$PROJECT_DIR/output/genomes/$TARGET_GENOME/annotations"
OUT_PREFIX="${TARGET_GENOME}-${TARGET_GENOME}-${SOURCE_ANNOT}"

# Check if output already exists
if [[ -f "$OUT_DIR/${OUT_PREFIX}.gtf.gz" ]] && [[ -f "$OUT_DIR/${OUT_PREFIX}.gff3.gz" ]]; then
    log "Output files already exist. Skipping."
    exit 0
fi

# Create output directory
mkdir -p "$OUT_DIR"
log "Starting liftoff for mapping $SOURCE_GENOME annotations to $TARGET_GENOME"
log "Working directory: $TEMP_DIR"

# Activate conda environment with liftoff
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate custom_genes
cd "$TEMP_DIR"

# Handle each input file
handle_input_file "$SOURCE_FA"
handle_input_file "$TARGET_FA" 
handle_input_file "$SOURCE_GTF"

# Set up liftoff command and add input files to command
log "Running liftoff..."
LIFTOFF_ARGS="-p $N_CORES -g $(basename "$SOURCE_GTF") -o ${OUT_PREFIX}.gff3 -dir $TEMP_DIR $(basename "$TARGET_FA" .gz) $(basename "$SOURCE_FA" .gz)"

# Check GTF and add -infer_genes if needed
gene_count=$(check_gene_feature "$(basename "$SOURCE_GTF")" | tail -1 )
if [ "$gene_count" -eq 0 ]; then
    log "Adding -infer_genes flag to liftoff command"
    LIFTOFF_ARGS="-infer_genes $LIFTOFF_ARGS"
fi

# Run liftoff with constructed arguments
log "Running liftoff with args: $LIFTOFF_ARGS"
liftoff $LIFTOFF_ARGS

# Convert GFF3 to GTF
log "Converting GFF3 to GTF..."
gffread "${OUT_PREFIX}.gff3" -F -T -o "${OUT_PREFIX}.gtf"

# Compress outputs
log "Compressing outputs..."
gzip "${OUT_PREFIX}.gff3" "${OUT_PREFIX}.gtf"

# Move results to final location
log "Moving results to output directory..."
rsync -av "${OUT_PREFIX}.gff3.gz" "${OUT_PREFIX}.gtf.gz" "$OUT_DIR/"

# Cleanup
log "Cleaning up temporary files..."
cd "$SCRATCH_DIR"
rm -rf "$TEMP_DIR"

log "Liftoff mapping complete"
conda deactivate
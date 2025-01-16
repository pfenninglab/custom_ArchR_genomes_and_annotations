#!/bin/bash
#SBATCH --partition=pool3-bigmem,pfen3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --error=logs/%x_%j.txt
#SBATCH --output=logs/%x_%j.txt

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

# Function to check if chromosome is a main chromosome
is_main_chromosome() {
    local chrom="$1"
    if [[ "$chrom" =~ ^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$ ]]; then
        return 0  # true
    else
        return 1  # false
    fi
}

# Function to filter GTF for main chromosomes
# assuming UCSC chromosome naming convention for human and mouse
filter_main_chromosomes() {
    local input_file="$1"
    local output_file="$2"
    
    log "Filtering GTF to keep only main chromosomes..."
    
    if [[ "$input_file" =~ \.gz$ ]]; then
        zcat "$input_file"
    else
        cat "$input_file"
    fi | awk '
        BEGIN { FS="\t"; OFS="\t" }
        /^#/ { print; next }
        $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ { print }
    ' > "$output_file"
    
    local orig_count=$(if [[ "$input_file" =~ \.gz$ ]]; then zcat "$input_file"; else cat "$input_file"; fi | grep -v "^#" | wc -l)
    local final_count=$(grep -v "^#" "$output_file" | wc -l)
    
    log "Filtered GTF from $orig_count to $final_count features"
    
    # Get stats on removed chromosomes
    log "Removed chromosome types:"
    if [[ "$input_file" =~ \.gz$ ]]; then
        zcat "$input_file"
    else
        cat "$input_file"
    fi | grep -v "^#" | cut -f1 | sort | uniq -c | 
    awk '! ($2 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/)' | head -n 10
}

# Check required tools
check_requirements() {
    local missing_tools=()
    for tool in samtools awk gffread zcat rsync; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing_tools+=("$tool")
        fi
    done
    
    if (( ${#missing_tools[@]} > 0 )); then
        log "ERROR: Missing required tools: ${missing_tools[*]}"
        exit 1
    fi
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
    
    if [ "$gene_count" -gt 0 ]; then
        log "Found $gene_count gene features in GTF"
    else
        log "No gene features found in first 1000 lines, will use -infer_genes"
    fi

    echo "$gene_count"
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

# Check required arguments and tools
if [[ -z "${SOURCE_GENOME:-}" ]] || [[ -z "${TARGET_GENOME:-}" ]] || [[ -z "${SOURCE_ANNOT:-}" ]]; then
    usage
fi

# Activate conda environment with liftoff
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate custom_genes
check_requirements

# Calculate resources
MEM_GB=$((${SLURM_MEM_PER_NODE:-64000} / 1024))
N_CORES=$((MEM_GB / 15))
log "Memory requested: ${MEM_GB}GB, using ${N_CORES} cores"

# Set up directories
PROJECT_DIR=${PROJECT_DIR:-$(pwd)}
SCRATCH_DIR=${SCRATCH_DIR:-/scratch/$USER}
OUT_DIR="$PROJECT_DIR/output/genomes/$TARGET_GENOME/annotations"
OUT_PREFIX="${TARGET_GENOME}-${SOURCE_GENOME}-${SOURCE_ANNOT}"
TEMP_DIR="$SCRATCH_DIR/${OUT_PREFIX}-${SLURM_JOB_ID}"

# Define paths
SOURCE_FA="$PROJECT_DIR/output/genomes/$SOURCE_GENOME/$SOURCE_GENOME.fa"
TARGET_FA="$PROJECT_DIR/output/genomes/$TARGET_GENOME/$TARGET_GENOME.fa"
SOURCE_GTF="$PROJECT_DIR/output/genomes/$SOURCE_GENOME/annotations/${SOURCE_ANNOT}.gtf"

# Check if output already exists
log "Checking if output files already exist..."
log "Looking for: $OUT_DIR/${OUT_PREFIX}.gtf.gz and $OUT_DIR/${OUT_PREFIX}.gff3.gz"

if [[ -f "$OUT_DIR/${OUT_PREFIX}.gtf.gz" && -f "$OUT_DIR/${OUT_PREFIX}.gff3.gz" ]]; then
    log "Output files already exist:"
    log "  - $OUT_DIR/${OUT_PREFIX}.gtf.gz"
    log "  - $OUT_DIR/${OUT_PREFIX}.gff3.gz"
    log "Skipping annotation mapping."
    exit 0
else
    log "Output files do not exist. Proceeding with annotation mapping."
fi

# Create directories
mkdir -p "$TEMP_DIR"/{input,filtered,output} "$OUT_DIR"
cd "$TEMP_DIR"

# Process input files
handle_input_file "$SOURCE_FA"
handle_input_file "$TARGET_FA" 
handle_input_file "$SOURCE_GTF"

# Get base filenames
SOURCE_FA_BASE=$(basename "$SOURCE_FA" .gz)
SOURCE_GTF_BASE=$(basename "$SOURCE_GTF" .gz)
TARGET_FA_BASE=$(basename "$TARGET_FA" .gz)

# Filter GTF for main chromosomes
FILTERED_GTF="$TEMP_DIR/filtered/${SOURCE_GTF_BASE%.gtf}.filtered.gtf"
filter_main_chromosomes "$SOURCE_GTF_BASE" "$FILTERED_GTF"

# Set up liftoff command
LIFTOFF_ARGS="-p ${N_CORES} -g ${FILTERED_GTF} -o ${OUT_PREFIX}.gff3 -dir ${TEMP_DIR}/output $TARGET_FA_BASE $SOURCE_FA_BASE"

# Check GTF and add -infer_genes if needed
gene_count=$(check_gene_feature "$FILTERED_GTF" | tail -1)
log "Gene count detected in GTF: $gene_count"

if [ "$gene_count" = "0" ]; then
    log "Adding -infer_genes flag to liftoff command"
    LIFTOFF_ARGS="-infer_genes $LIFTOFF_ARGS"
fi

# Run liftoff
log "Running liftoff with args: $LIFTOFF_ARGS"
liftoff $LIFTOFF_ARGS

# Convert and compress outputs
log "Converting GFF3 to GTF..."
gffread "${OUT_PREFIX}.gff3" -F -T -o "${OUT_PREFIX}.gtf"

log "Compressing outputs..."
gzip -f "${OUT_PREFIX}.gff3" "${OUT_PREFIX}.gtf"

# Move results to final location
log "Moving results to output directory..."
rsync -av --remove-source-files ${OUT_PREFIX}.*.gz $OUT_DIR/

# Clean up
rm -rf "$TEMP_DIR"

log "Liftoff mapping complete"
conda deactivate

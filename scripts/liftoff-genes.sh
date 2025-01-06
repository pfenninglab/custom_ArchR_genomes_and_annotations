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

# Function to get chromosomes from FASTA and filter GTF
filter_gtf_by_fasta_chroms() {
    local fasta_file="$1"
    local gtf_file="$2"
    local output_gtf="$3"
    local trackname="$4"
       
    log "Indexing FASTA file to get chromosomes..."
    samtools faidx "$fasta_file"
    
    # Create a temporary directory for processing
    local tmp_dir="tmp_filter_$$"
    mkdir -p "$tmp_dir"
    
    # Get chromosomes from FASTA index into a file
    cut -f1 "${fasta_file}.fai" | sort > "$tmp_dir/fasta_chroms.txt"
    local chrom_count=$(wc -l < "$tmp_dir/fasta_chroms.txt")
    log "Found $chrom_count chromosomes in FASTA"
    
    # Get unique chromosomes from GTF into a file
    if [[ "$gtf_file" =~ \.gz$ ]]; then
        zcat "$gtf_file" | grep -v "^#" | cut -f1 | sort -u > "$tmp_dir/gtf_chroms.txt"
    else
        grep -v "^#" "$gtf_file" | cut -f1 | sort -u > "$tmp_dir/gtf_chroms.txt"
    fi
    local gtf_chrom_count=$(wc -l < "$tmp_dir/gtf_chroms.txt")
    log "Found $gtf_chrom_count chromosomes in GTF"
    
    # Find common chromosomes
    comm -12 "$tmp_dir/fasta_chroms.txt" "$tmp_dir/gtf_chroms.txt" > "$tmp_dir/common_chroms.txt"
    local common_count=$(wc -l < "$tmp_dir/common_chroms.txt")
    log "Found $common_count chromosomes in common between FASTA and GTF"

    # Create a filtered GTF with only common chromosomes
    if [[ "$gtf_file" =~ \.gz$ ]]; then
        log "Filtering compressed GTF..."
        zcat "$gtf_file" | awk -v chroms="$tmp_dir/common_chroms.txt" '
            BEGIN { while((getline chr < chroms) > 0) valid_chroms[chr]=1 }
            /^#/ { print; next }
            valid_chroms[$1] { print }
        ' > "$output_gtf"
    else
        log "Filtering uncompressed GTF..."
        awk -v chroms="$tmp_dir/common_chroms.txt" '
            BEGIN { while((getline chr < chroms) > 0) valid_chroms[chr]=1 }
            /^#/ { print; next }
            valid_chroms[$1] { print }
        ' "$gtf_file" > "$output_gtf"
    fi

    # Log feature counts
    local orig_count
    local final_count
    if [[ "$gtf_file" =~ \.gz$ ]]; then
        orig_count=$(zcat "$gtf_file" | grep -v "^#" | wc -l)
    else
        orig_count=$(grep -v "^#" "$gtf_file" | wc -l)
    fi
    final_count=$(grep -v "^#" "$output_gtf" | wc -l)
    
    log "Original GTF had $orig_count features"
    log "Filtered GTF has $final_count features"
    
    # Validate output
    if [[ ! -s "$output_gtf" ]]; then
        log "ERROR: Filtered GTF is empty after chromosome filtering"
        rm -rf "$tmp_dir"
        exit 1
    fi
    
    # Print some example chromosomes that were filtered out
    log "Example chromosomes present in GTF but not in FASTA:"
    comm -23 "$tmp_dir/gtf_chroms.txt" "$tmp_dir/fasta_chroms.txt" | head -n 5 | while read -r chrom; do
        log "  - $chrom"
    done

    # Cleanup
    rm -rf "$tmp_dir"
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

# Activate conda environment with liftoff
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate custom_genes
check_requirements

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
OUT_PREFIX="${TARGET_GENOME}_liftoff_${SOURCE_ANNOT}"

# Check if output already exists
if [[ -f "$OUT_DIR/${OUT_PREFIX}.gtf.gz" ]] && [[ -f "$OUT_DIR/${OUT_PREFIX}.gff3.gz" ]]; then
    log "Output files already exist. Skipping."
    exit 0
fi

# Create output directory
mkdir -p "$OUT_DIR"
log "Starting liftoff for mapping $SOURCE_GENOME annotations to $TARGET_GENOME"
log "Working directory: $TEMP_DIR"

cd "$TEMP_DIR"

# Handle each input file
handle_input_file "$SOURCE_FA"
handle_input_file "$TARGET_FA" 
handle_input_file "$SOURCE_GTF"

# Get base filenames without .gz extension
SOURCE_FA_BASE=$(basename "$SOURCE_FA" .gz)
SOURCE_GTF_BASE=$(basename "$SOURCE_GTF" .gz)
TARGET_FA_BASE=$(basename "$TARGET_FA" .gz)

# Create subdirectories for better organization
mkdir -p "$TEMP_DIR"/{input,filtered,output}

# Filter GTF to match FASTA chromosomes
log "Filtering GTF to match FASTA chromosomes..."
FILTERED_GTF="$TEMP_DIR/filtered/${SOURCE_GTF_BASE%.gtf}.filtered.gtf"
filter_gtf_by_fasta_chroms "$SOURCE_FA_BASE" "$SOURCE_GTF_BASE" "$FILTERED_GTF" "$SOURCE_ANNOT"

# Validate filtered GTF
if [[ ! -s "$FILTERED_GTF" ]]; then
    log "ERROR: Filtered GTF is empty. Check that FASTA and GTF have matching chromosomes."
    exit 1
fi

# Set up liftoff command
log "Running liftoff..."
LIFTOFF_ARGS="-p $N_CORES -g $FILTERED_GTF -o ${OUT_PREFIX}.gff3 -dir $TEMP_DIR/output $TARGET_FA_BASE $SOURCE_FA_BASE"

# Check GTF and add -infer_genes if needed 
gene_count=$(check_gene_feature "$FILTERED_GTF" | tail -n 1)
if [ "$gene_count" -eq 0 ]; then
    log "Adding -infer_genes flag to liftoff command"
    LIFTOFF_ARGS="-infer_genes $LIFTOFF_ARGS"
fi

# Run liftoff with constructed arguments
log "Running liftoff with args: $LIFTOFF_ARGS"
liftoff $LIFTOFF_ARGS

# Verify liftoff output exists
if [[ ! -f "${OUT_PREFIX}.gff3" ]]; then
    log "ERROR: Liftoff failed to produce output GFF3 file"
    exit 1
fi

# Convert GFF3 to GTF
log "Converting GFF3 to GTF..."
if ! gffread "${OUT_PREFIX}.gff3" -F -T -o "${OUT_PREFIX}.gtf"; then
    log "ERROR: Failed to convert GFF3 to GTF"
    exit 1
fi

# Verify GTF conversion
if [[ ! -s "${OUT_PREFIX}.gtf" ]]; then
    log "ERROR: GTF conversion produced empty file"
    exit 1
fi

# Compress outputs
log "Compressing outputs..."
for file in "${OUT_PREFIX}.gff3" "${OUT_PREFIX}.gtf"; do
    if ! gzip -f "$file"; then
        log "ERROR: Failed to compress $file"
        exit 1
    fi
done

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Move results to final location
log "Moving results to output directory..."
for file in "${OUT_PREFIX}.gff3.gz" "${OUT_PREFIX}.gtf.gz"; do
    if ! rsync -av "$file" "$OUT_DIR/"; then
        log "ERROR: Failed to copy $file to output directory"
        exit 1
    fi
done

# Verify files were copied successfully
for file in "gff3.gz" "gtf.gz"; do
    if [[ ! -f "$OUT_DIR/${OUT_PREFIX}.$file" ]]; then
        log "ERROR: Failed to find copied file: $OUT_DIR/${OUT_PREFIX}.$file"
        exit 1
    fi
done

log "Liftoff mapping complete"
conda deactivate
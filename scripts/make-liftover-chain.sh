#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --time=3-00:00:00
#SBATCH --export=ALL
#SBATCH --mem=24G
#SBATCH --error=logs/%x_%j.txt
#SBATCH --output=logs/%x_%j.txt

set -euo pipefail

# Function to print usage
usage() {
    echo "Usage: $0 -k KENT_DIR -t TARGET -q QUERY -c1 TARGET_CLADE -c2 QUERY_CLADE [-o OUTPUT_DIR] [-w SCRATCH_DIR]" 1>&2
    echo ""
    echo "Required arguments:"
    echo "  -k KENT_DIR    Path to Kent source repo directory"
    echo "  -t TARGET      Target genome FASTA file (.fa or .fa.gz)"
    echo "  -q QUERY       Query genome FASTA file (.fa or .fa.gz)"
    echo "  -c1 CLADE1     Target clade (one of: mammal, primate, vertebrate)"
    echo "  -c2 CLADE2     Query clade (one of: mammal, primate, vertebrate)"
    echo ""
    echo "Optional arguments:"
    echo "  -o OUTPUT_DIR  Output directory (default: \$PWD)"
    echo "  -w SCRATCH_DIR Scratch directory (default: /scratch/\$USER)"
    echo "  -h            Show this help message"
    exit 1
}

# Function to log messages with timestamps 
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to get base name without extensions
get_base_name() {
    local filename="$1"
    basename="${filename##*/}"
    basename="${basename%.gz}"
    basename="${basename%.fa}"
    basename="${basename%.fasta}"
    echo "$basename"
}

# Function to check required software
check_requirements() {
    local required_tools=(
        "faToTwoBit"
        "twoBitInfo"
    )
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            log "ERROR: Required tool '$tool' not found in PATH"
            exit 1
        fi
    done
    
    # Check for pairLastz.sh in Kent directory
    local pairLastz="${KENT_DIR}/src/hg/utils/automation/pairLastz.sh"
    if [[ ! -f "$pairLastz" ]]; then
        log "ERROR: pairLastz.sh not found or not executable at $pairLastz"
        exit 1
    fi
}

# Function to process FASTA to 2bit
process_to_2bit() {
    local input_file="$1"
    local output_dir="$2"
    local genome_name=$(get_base_name "$input_file")
    local genome_dir=$(dirname "$input_file")
    local existing_2bit="${genome_dir}/${genome_name}.2bit"
    
    mkdir -p "$output_dir"
    
    if [[ -f "$existing_2bit" ]]; then
        log "Found existing 2bit file for ${genome_name}, copying to scratch..."
        rsync -P "$existing_2bit" "${output_dir}/${genome_name}.2bit"
    else
        if [[ "$input_file" =~ \.gz$ ]]; then
            log "Processing gzipped file ${genome_name}..."
            gunzip -c "$input_file" > "${output_dir}/${genome_name}.fa"
        else
            log "Processing uncompressed file ${genome_name}..."
            rsync -Paq "$input_file" "${output_dir}/${genome_name}.fa"
        fi
        
        log "Converting ${genome_name} to 2bit format..."
        faToTwoBit "${output_dir}/${genome_name}.fa" "${output_dir}/${genome_name}.2bit"
        
        log "Copying 2bit file back to genome directory..."
        rsync -P "${output_dir}/${genome_name}.2bit" "$existing_2bit"
        
        rm "${output_dir}/${genome_name}.fa"
    fi
    
    # Generate chromosome sizes
    twoBitInfo "${output_dir}/${genome_name}.2bit" "${output_dir}/${genome_name}.chrom.sizes"
    
    echo "${output_dir}/${genome_name}.2bit"
}

# Function to validate clade
validate_clade() {
    local clade="$1"
    local valid_clades=("mammal" "primate" "vertebrate")
    
    for valid_clade in "${valid_clades[@]}"; do
        if [[ "$clade" == "$valid_clade" ]]; then
            return 0
        fi
    done
    
    log "ERROR: Invalid clade '$clade'. Must be one of: ${valid_clades[*]}"
    return 1
}

# Parse command line arguments
while getopts "k:t:q:c1:c2:o:w:h" opt; do
    case $opt in
        k) KENT_DIR="$OPTARG" ;;
        t) TARGET="$OPTARG" ;;
        q) QUERY="$OPTARG" ;;
        c1) TARGET_CLADE="$OPTARG" ;;
        c2) QUERY_CLADE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        w) SCRATCH_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [[ -z "${KENT_DIR:-}" || -z "${TARGET:-}" || -z "${QUERY:-}" || \
      -z "${TARGET_CLADE:-}" || -z "${QUERY_CLADE:-}" ]]; then
    usage
fi

# Ensure Kent directory exists and contains required script
if [[ ! -d "$KENT_DIR" ]]; then
    log "ERROR: Kent directory not found at $KENT_DIR"
    exit 1
fi

# Validate clades
validate_clade "$TARGET_CLADE" || exit 1
validate_clade "$QUERY_CLADE" || exit 1

# Set default directories if not specified
OUTPUT_DIR=${OUTPUT_DIR:-"$PWD"}
SCRATCH_DIR=${SCRATCH_DIR:-"/scratch/$USER"}

# Check requirements
conda activate custom_genes
check_requirements

# Check input files exist and are readable
for file in "$TARGET" "$QUERY"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: File not found at $file"
        exit 1
    fi
    if [[ ! -r "$file" ]]; then
        log "ERROR: File exists but is not readable at $file"
        exit 1
    fi
done

# Extract genome names
TARGET_NAME=$(get_base_name "$TARGET")
QUERY_NAME=$(get_base_name "$QUERY")

# Set up working directories
TEMP_DIR="${SCRATCH_DIR}/lastz_${TARGET_NAME}_${QUERY_NAME}_${SLURM_JOB_ID}"
mkdir -p "$TEMP_DIR"/{input,temp,output} "$OUTPUT_DIR"
log "Working in temporary directory: $TEMP_DIR"

# Process genomes to 2bit format
TARGET_2BIT=$(process_to_2bit "$TARGET" "$TEMP_DIR/input"| tail -n 1)
QUERY_2BIT=$(process_to_2bit "$QUERY" "$TEMP_DIR/input"| tail -n 1)

# Run pairLastz alignment
log "Running LASTZ alignment..."
cd "$TEMP_DIR"

# Use the validated clades from command line arguments
log "Using target clade: $TARGET_CLADE, query clade: $QUERY_CLADE"

${KENT_DIR}/src/hg/utils/automation/pairLastz.sh \
"$FINAL_TARGET" "$FINAL_QUERY" "$TARGET_CLADE" "$QUERY_CLADE"

# Copy results back to output directory
log "Copying results to output directory..."
rsync -av "$TEMP_DIR/" "$OUTPUT_DIR/"
rsync -av "$TEMP_DIR/lastz.log" "$OUTPUT_DIR/"

# Cleanup
log "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

log "LASTZ alignment complete. Results are in $OUTPUT_DIR"
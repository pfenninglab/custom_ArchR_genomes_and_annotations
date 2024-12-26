#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --time 1-00:00:00
#SBATCH --export=ALL
#SBATCH --mem=8G
#SBATCH --error=logs/%x_%j.txt
#SBATCH --output=logs/%x_%j.txt

set -euo pipefail

# Function to print usage
usage() {
    echo "Usage: $0 -g GENOME -f FASTA_URL [-p PROJECT_DIR] [-s SCRATCH_DIR] [-n GTF_NAME -t GTF_URL] [-b BLACKLIST_URL]" 1>&2
    echo "Downloads and prepares genome files for analysis pipeline"
    echo ""
    echo "Required arguments:"
    echo "  -g GENOME        Genome name/version (e.g. hg38)"
    echo "  -f FASTA_URL     URL to genome FASTA file"
    echo ""
    echo "Optional arguments:"
    echo "  -p PROJECT_DIR   Project directory (default: current directory)"
    echo "  -s SCRATCH_DIR   Scratch directory for computations"
    echo "  -n GTF_NAME      Name for GTF annotation (e.g. gencode.v39.basic)"
    echo "  -t GTF_URL       URL to gene annotation GTF file"
    echo "  -b BLACKLIST_URL URL to blacklist file"
    echo "  -h               Show this help message"
    exit 1
}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to extract GTF name from URL
get_name() {
    local url="$1"
    local basename=$(basename "$url")
    # Strip .gz extension if present
    basename=${basename%.gz}
    # Strip .gtf extension
    basename=${basename%.gtf}
    echo "$basename"
}

# Function to check if file exists in either scratch or project directory
check_file() {
    local basename="$1"
    local dir="$2"
    local scratch_path="$SCRATCH_DIR/$dir/$basename"
    local project_path="$PROJECT_DIR/output/genomes/$dir/$basename"
    
    [[ -f "$scratch_path" ]] || [[ -f "$project_path" ]]
}

# Function to sync file from scratch to project directory
sync_file() {
    local source_file="$1"
    local relative_path="${source_file#$SCRATCH_DIR/}"
    local target_dir="$PROJECT_DIR/output/genomes/$(dirname "$relative_path")"
    
    if [[ "$SCRATCH_DIR" != "$PROJECT_DIR" ]]; then
        log "Syncing file to project directory: $(basename "$source_file")"
        mkdir -p "$target_dir"
        rsync -Pav "$source_file" "$target_dir/"
    fi
}

# Function to download file using wget or rsync depending on URL
download_file() {
    local url="$1"
    local output_file="$2"
    local dir="$3"
    local output_dir=$(dirname "$output_file")
    
    # Check if compressed file exists
    if check_file "$(basename "$output_file")" "$dir"; then
        log "File already exists: $(basename "$output_file")"
        return 0
    fi
    
    mkdir -p "$output_dir"
    
    if [[ "$url" == *"hgdownload.soe.ucsc.edu"* ]]; then
        # Convert UCSC http URL to rsync URL
        local rsync_url=${url/https:/rsync:}
        log "Using rsync to download: $rsync_url"
        rsync -avP "$rsync_url" "$output_file"
    else
        log "Downloading: $url"
        wget --tries=0 --no-check-certificate -O "$output_file" "$url"
    fi
    
    # Sync file immediately after download if using scratch
    sync_file "$output_file"
}

# Parse command line arguments
while getopts "g:f:n:t:b:p:s:h" opt; do
    case $opt in
        g) GENOME="$OPTARG" ;;
        f) FASTA_URL="$OPTARG" ;;
        n) GTF_NAME="$OPTARG" ;;
        t) GTF_URL="$OPTARG" ;;
        b) BLACKLIST_URL="$OPTARG" ;;
        p) PROJECT_DIR="$OPTARG" ;;
        s) SCRATCH_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# If GENOME is not provided, extract name from the URL
if [[ -z "${GENOME:-}" ]] || [[ -z "${FASTA_URL:-}" ]]; then
    if [[ -z "${GENOME:-}" ]]; then
        GENOME=$(get_name "$FASTA_URL")
        log "Auto-detected GENOME name from URL: $FASTA_URL"
    fi
fi

# If GTF URL is provided but no name, extract name from URL
if [[ -n "${GTF_URL:-}" ]]; then
    if [[ -z "${GTF_NAME:-}" ]]; then
        GTF_NAME=$(get_name "$GTF_URL")
        log "Auto-detected GTF name from URL: $GTF_NAME"
    fi
fi

# Set default directories if not specified
PROJECT_DIR=${PROJECT_DIR:-$(pwd)}
SCRATCH_DIR=${SCRATCH_DIR:-$PROJECT_DIR}

# Create directory structure
GENOME_DIR="$SCRATCH_DIR/$GENOME"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$GENOME_DIR" "$LOG_DIR"

log "Starting genome download for $GENOME"
log "Project directory: $PROJECT_DIR"
log "Scratch directory: $SCRATCH_DIR"

# Download FASTA
FASTA_FILE="$GENOME_DIR/$GENOME.fa.gz"
download_file "$FASTA_URL" "$FASTA_FILE" "$GENOME"

# Download GTF if URL provided
if [[ -n "${GTF_URL:-}" ]]; then
    GTF_DIR="$GENOME_DIR/annotations"
    mkdir -p "$GTF_DIR"
    GTF_FILE="$GTF_DIR/${GTF_NAME}.gtf.gz"
    download_file "$GTF_URL" "$GTF_FILE" "annotations"
fi

# Download blacklist if URL provided
if [[ -n "${BLACKLIST_URL:-}" ]]; then
    BLACKLIST_FILE="$GENOME_DIR/${GENOME}_blacklist.bed.gz"
    download_file "$BLACKLIST_URL" "$BLACKLIST_FILE" "$GENOME"
fi

log "Genome download and preparation complete"
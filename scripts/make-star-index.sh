#!/bin/bash
#!/bin/bash
#SBATCH --partition=pool3-bigmem,pfen3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --error=logs/%x_%j.txt
#SBATCH --output=logs/%x_%j.txt

#
# Build STAR genome index with proper scratch directory handling and error checking
# Supports both gzipped and uncompressed input files
#

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Function to print usage
usage() {
    echo "Usage: $0 -f FASTA -t GTF -o OUTPUT_DIR [-p THREADS] [-w SCRATCH_DIR] [-m MEM]" 1>&2
    echo ""
    echo "Required arguments:"
    echo "  -f FASTA       Path to genome FASTA file (can be gzipped)"
    echo "  -t GTF         Path to gene annotations GTF file (can be gzipped)"
    echo "  -o OUTPUT_DIR  Output directory for STAR index"
    echo ""
    echo "Optional arguments:"
    echo "  -p THREADS     Number of threads (default: derived from SLURM_CPUS_PER_TASK)"
    echo "  -w SCRATCH_DIR Scratch directory (default: /scratch/\$USER)"
    echo "  -m MEM        Memory in GB (default: 90% of SLURM_MEM_PER_NODE)"
    echo "  -h            Show this help message"
    exit 1
}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to handle input files
handle_input_file() {
    local file="$1"
    local base=$(basename "$file" .gz)
    
    if [[ "$file" =~ \.gz$ ]]; then
        log "Using compressed file: $file"
        rsync -Paq "$file" "./${base}.gz"
        gunzip -f "./${base}.gz"
    else
        log "Using uncompressed file: $file"
        rsync -Paq "$file" "./${base}"
    fi
}

# Parse command line arguments
while getopts "g:f:t:o:p:w:m:h" opt; do
    case $opt in
        f) FASTA="$OPTARG" ;;
        t) GTF="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        w) SCRATCH_DIR="$OPTARG" ;;
        m) MEM="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate custom_genes

# Check required arguments
if [[ -z "${FASTA:-}" || -z "${GTF:-}" || -z "${OUTPUT_DIR:-}" ]]; then
    usage
fi

# Set default values for optional arguments
MEM_GB=$((${SLURM_MEM_PER_NODE:-64000} / 1024))
THREADS=$((MEM_GB / 15))

SCRATCH_DIR=${SCRATCH_DIR:-/scratch/$USER}
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    # Convert SLURM memory from MB to GB and take 90%
    MEM=${MEM:-$(( SLURM_MEM_PER_NODE * 90 / 1024 / 100 ))}
else
    MEM=${MEM:-32}
fi

# Validate input files exist (either compressed or uncompressed)
for file in "$FASTA" "$GTF"; do
    if [[ ! -f "$file" && ! -f "${file}.gz" ]]; then
        log "ERROR: Neither $file nor ${file}.gz found"
        exit 1
    fi
done

# Set up directories
TEMP_DIR="${SCRATCH_DIR}/star_${SLURM_JOB_ID}"
TEMP_INPUT="${TEMP_DIR}/input"
TEMP_OUTPUT="${TEMP_DIR}/output"

log "Creating temporary directories..."
mkdir -p "$TEMP_INPUT" "$TEMP_OUTPUT"

# Function to clean up scratch directory
cleanup() {
    local exit_code=$?
    log "Cleaning up scratch directory..."
    if [ -d "$TEMP_DIR" ]; then
        rm -rf "$TEMP_DIR"
    fi
    exit $exit_code
}

# Set up cleanup trap
trap cleanup EXIT
# Change to input directory for file handling
cd "$TEMP_INPUT"

# Handle input files
log "Processing input files..."
handle_input_file "$FASTA"
handle_input_file "$GTF"

# Get base filenames without .gz extension
FASTA_BASE=$(basename "$FASTA" .gz)
GTF_BASE=$(basename "$GTF" .gz)

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Validate read length parameter (for sjdbOverhang)
READ_LENGTH=${READ_LENGTH:-100}  # Default to 100bp reads
OVERHANG=$(( READ_LENGTH - 1 ))

# Build STAR index with updated parameters
log "Building STAR index with ${THREADS} threads and ${MEM}GB memory..."
STAR --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$TEMP_OUTPUT" \
    --genomeFastaFiles "$TEMP_INPUT/$FASTA_BASE" \
    --sjdbGTFfile "$TEMP_INPUT/$GTF_BASE" \
    --sjdbOverhang "$OVERHANG" \
    --limitGenomeGenerateRAM $(( MEM * 1024**3 * 9 / 10 ))

# Check if index was created successfully
if [[ ! -f "${TEMP_OUTPUT}/Genome" ]]; then
    log "ERROR: STAR index generation failed"
    exit 1
fi

# Copy results back to output directory
log "Copying STAR index to output directory..."
rsync -Paq --remove-source-files ${TEMP_OUTPUT}/ ${OUTPUT_DIR}/

rm -r $TEMP_DIR
log ""

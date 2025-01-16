#!/bin/bash
#
# Execute STAR index generation for all valid source-target genome combinations
#

PROJECT_DIR=${1:-$HOME/repos/custom_ArchR_genomes_and_annotations}
SCRATCH_DIR=${2:-/scratch/$USER}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Process source genomes config
log "Reading source genomes config..."
config_file="$PROJECT_DIR/config/source_genomes.tsv"
num_sources=$(wc -l < "$config_file")

# Get list of target genomes
log "Reading target genomes config..."
targets_file="$PROJECT_DIR/config/target_genomes.tsv"
num_targets=$(wc -l < "$targets_file")

# Create logs directory
mkdir -p "$PROJECT_DIR/logs"

# Process human genomes to non-rodents (lines 2-5 of source config)
i=2    # Parse source genome info
row=$(awk -v line=$i 'NR==line' "$config_file")
SOURCE_GENOME=$(echo "$row" | awk -F'\t' '{print $2}')
SOURCE_ANNOT=$(echo "$row" | awk -F'\t' '{print $5}' | tr -d '\r\n')

# Loop through non-rodent target genomes (lines 2-8 of target config)
for j in $(seq 2 10); do
    # Parse target genome info
    row=$(awk -v line=$j 'NR==line' "$targets_file")
    TARGET_GENOME=$(echo "$row" | awk -F'\t' '{print $2}' | tr -d '\r\n')
    
    # Construct paths
    OUT_PREFIX="${TARGET_GENOME}-${SOURCE_GENOME}-${SOURCE_ANNOT}"
    GTF="$PROJECT_DIR/output/genomes/$TARGET_GENOME/annotations/${OUT_PREFIX}.gtf.gz"
    FASTA="$PROJECT_DIR/output/genomes/$TARGET_GENOME/$TARGET_GENOME.fa.gz"
    OUTPUT_DIR="$PROJECT_DIR/output/genomes/$TARGET_GENOME/star_index/${OUT_PREFIX}"

    # Check if files exist and STAR index doesn't
    if [[ -f "$GTF" && -f "$FASTA" && ! -f "$OUTPUT_DIR/Genome" ]]; then
        log "Submitting STAR index job for: $OUT_PREFIX"
        
        # Submit SLURM job
        sbatch --job-name="star_${OUT_PREFIX}" \
            --partition pool3-bigmem --mem 64G \
            --time 24:00:00 \
            "$PROJECT_DIR/scripts/make-star-index.sh" \
            -f "$FASTA" \
            -t "$GTF" \
            -o "$OUTPUT_DIR"
    else
        if [[ ! -f "$GTF" ]]; then
            log "WARNING: Missing lifted GTF: $GTF"
        fi
        if [[ ! -f "$FASTA" ]]; then
            log "WARNING: Missing genome FASTA: $FASTA"
        fi
        if [[ -d "$OUTPUT_DIR" ]]; then
            log "STAR index already exists: $OUTPUT_DIR"
        fi
    fi
done



# Process mouse genomes to rodents (lines 6-8 of source config)
i=6    # Parse source genome info
row=$(awk -v line=$i 'NR==line' "$config_file")
SOURCE_GENOME=$(echo "$row" | awk -F'\t' '{print $2}')
SOURCE_ANNOT=$(echo "$row" | awk -F'\t' '{print $5}' | tr -d '\r\n')

# Loop through rodent target genomes (lines 9-10 of target config)
for j in $(seq 11 12); do
    # Parse target genome info
    row=$(awk -v line=$j 'NR==line' "$targets_file")
    TARGET_GENOME=$(echo "$row" | awk -F'\t' '{print $2}' | tr -d '\r\n')
    
    # Construct paths
    OUT_PREFIX="${TARGET_GENOME}-${SOURCE_GENOME}-${SOURCE_ANNOT}"
    GTF="$PROJECT_DIR/output/genomes/$TARGET_GENOME/annotations/${OUT_PREFIX}.gtf.gz"
    FASTA="$PROJECT_DIR/output/genomes/$TARGET_GENOME/$TARGET_GENOME.fa.gz"
    OUTPUT_DIR="$PROJECT_DIR/output/genomes/$TARGET_GENOME/star_index/${OUT_PREFIX}"

    # Check if files exist and STAR index doesn't
    if [[ -f "$GTF" && -f "$FASTA" && ! -f "$OUTPUT_DIR/Genome" ]]; then
        log "Submitting STAR index job for: $OUT_PREFIX"
        
        sbatch --job-name="star_${OUT_PREFIX}" \
            --partition pool3-bigmem --mem 64G \
            --time 24:00:00 \
            "$PROJECT_DIR/scripts/make-star-index.sh" \
            -f "$FASTA" \
            -t "$GTF" \
            -o "$OUTPUT_DIR"
    else
        if [[ ! -f "$GTF" ]]; then
            log "WARNING: Missing lifted GTF: $GTF"
        fi
        if [[ ! -f "$FASTA" ]]; then
            log "WARNING: Missing genome FASTA: $FASTA"
        fi
        if [[ -d "$OUTPUT_DIR" ]]; then
            log "STAR index already exists: $OUTPUT_DIR"
        fi
    fi
done

log "All STAR index jobs submitted"

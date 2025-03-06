#!/bin/bash
# File: inst/scripts/run-all-archr-genomes.sh

set -euo pipefail

PROJECT_DIR=${1:-$(pwd)}
SCRATCH_DIR=${2:-/scratch/$USER}

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log "Starting ArchR genome generation for all target genomes..."

# Read target genomes config
while IFS=$'\t' read -r common genome species url liftover_from_hg38 liftover_from_mm10 liftover_to_h38 liftover_to_mm10 bsgenome extra; do
    # Skip header
    [[ "$common" == "common_name" ]] && continue
    
    # Check if genome directory exists
    if [[ ! -d "$PROJECT_DIR/genomes/$genome" ]]; then
        log "Skipping $genome - directory not found"
        continue
    fi
    
    # Process human annotations
    if [[ -n "$liftover_from_hg38" ]]; then
        for annot in "gencode.v47.basic" "gencode.v47.comp" "gencode.v44.basic"; do
            if [[ -f "$PROJECT_DIR/genomes/$genome/annotations/${genome}_liftoff_hg38_${annot}.gtf.gz" ]]; then
                log "Submitting job for $genome with hg38 $annot annotation"
                sbatch --job-name="archr_${genome}_${annot}" \
                    "$PROJECT_DIR/inst/scripts/create-archr-genome.sh" \
                    -g "$genome" \
                    -s "hg38" \
                    -a "$annot" \
                    -p "$PROJECT_DIR" \
                    -w "$SCRATCH_DIR"
            fi
        done
    fi
    
    # Process mouse annotations
    if [[ -n "$liftover_from_mm10" ]]; then
        for annot in "gencode.vM25.basic" "gencode.vM25.comp"; do
            if [[ -f "$PROJECT_DIR/genomes/$genome/annotations/${genome}_liftoff_mm10_${annot}.gtf.gz" ]]; then
                log "Submitting job for $genome with mm10 $annot annotation"
                sbatch --job-name="archr_${genome}_${annot}" \
                    "$PROJECT_DIR/inst/scripts/create-archr-genome.sh" \
                    -g "$genome" \
                    -s "mm10" \
                    -a "$annot" \
                    -p "$PROJECT_DIR" \
                    -w "$SCRATCH_DIR"
            fi
        done
    fi
done < "$PROJECT_DIR/inst/config/target_genomes.tsv"

log "All ArchR genome generation jobs submitted"

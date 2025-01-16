#!/bin/bash
#
# Execute Chromap index builds from config file
#

PROJECT_DIR=$HOME/repos/custom_ArchR_genomes_and_annotations
SCRATCH_DIR=/scratch/$USER
GENOME_DIR=${PROJECT_DIR}/output/genomes

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Process genomes from config file
conda activate custom_genes

log "Processing genomes..."
config_file="config/target_genomes.tsv"
num_genomes=$(wc -l < "$config_file")

# Process each genome starting from line 2 (skip header)
for i in $(seq 2 12); do
    # Parse TSV line using awk
    row=$(awk -v line="$i" 'NR==line' "$config_file")
    genome=$(echo "$row" | awk -F'\t' '{print $2}')
    fasta_path="${GENOME_DIR}/${genome}/${genome}.fa.gz"
    output_dir="${GENOME_DIR}/${genome}/chromap_index"
    
    # Skip if index already exists
    if [[ -f "${output_dir}/index" ]]; then
        log "Chromap index already exists for ${genome}, skipping..."
        continue
    fi
    
    # Check if FASTA exists
    if [[ ! -f "$fasta_path" ]]; then
        log "WARNING: FASTA file not found for ${genome} at ${fasta_path}, skipping..."
        continue
    fi
    
    log "Submitting job for genome: $genome"
    sbatch --job-name="chromap_${genome}" scripts/make-chromap-index.sh \
        -g "$genome" \
        -f "$fasta_path" \
        -o "$output_dir" \
        -w "$SCRATCH_DIR"
done


# do this for the human and mouse genomes
config_file="config/source_genomes.tsv"
# Process each genome starting from line 2 (skip header)
for i in 2 6; do
    # Parse TSV line using awk
    row=$(awk -v line="$i" 'NR==line' "$config_file")
    genome=$(echo "$row" | awk -F'\t' '{print $2}')
    fasta_path="${GENOME_DIR}/${genome}/${genome}.fa.gz"
    output_dir="${GENOME_DIR}/${genome}/chromap_index"
    
    # Skip if index already exists
    if [[ -f "${output_dir}/index" ]]; then
        log "Chromap index already exists for ${genome}, skipping..."
        continue
    fi
    
    # Check if FASTA exists
    if [[ ! -f "$fasta_path" ]]; then
        log "WARNING: FASTA file not found for ${genome} at ${fasta_path}, skipping..."
        continue
    fi
    
    log "Submitting job for genome: $genome"
    sbatch --job-name="chromap_${genome}" scripts/make-chromap-index.sh \
        -g "$genome" \
        -f "$fasta_path" \
        -o "$output_dir" \
        -w "$SCRATCH_DIR"
done
log "All Chromap index build jobs submitted"




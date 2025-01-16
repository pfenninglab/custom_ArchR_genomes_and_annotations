#!/bin/bash
#
# Execute genome downloads for source and target genomes
#

PROJECT_DIR=$HOME/repos/custom_ArchR_genomes_and_annotations
SCRATCH_DIR=/scratch/$USER

# Process source genomes
echo "Processing source genomes..."
config_file="config/source_genomes.tsv"
num_sources=$(wc -l < $config_file)
conda activate custom_genes

for i in $(seq 2 $num_sources); do
    # Parse TSV line using awk
    row=$(awk -v line=$i 'NR==line' $config_file)
    genome=$(echo "$row" | awk -F'\t' '{print $2}')
    fasta_url=$(echo "$row" | awk -F'\t' '{print $4}')
    gtf_name=$(echo "$row" | awk -F'\t' '{print $5}')
    gtf_url=$(echo "$row" | awk -F'\t' '{print $6}')
    # last column trim the newline character
    blacklist_url=$(echo "$row" | awk -F'\t' '{print $7}'| tr -d '\r\n')

    echo "Submitting job for source genome: $genome"
    sbatch --job-name="dl_${genome}" scripts/download-genome.sh \
        -g "$genome" \
        -f "$fasta_url" \
        -n "$gtf_name" \
        -t "$gtf_url" \
        -b "$blacklist_url" \
        -p "$PROJECT_DIR" \
        -s "$SCRATCH_DIR"
done

# Process target genomes 
echo "Processing target genomes..."
config_file="config/target_genomes.tsv"
num_targets=$(wc -l < $config_file)
for i in $(seq 2 12); do
    # Parse TSV line using awk
    row=$(awk -v line=$i 'NR==line' $config_file)
    genome=$(echo "$row" | awk -F'\t' '{print $2}')
    fasta_url=$(echo "$row" | awk -F'\t' '{print $4}')
    
    # Extract liftover chain URLs and combine them
    liftover_from_hg38=$(echo "$row" | awk -F'\t' '{print $5}')
    liftover_from_mm10=$(echo "$row" | awk -F'\t' '{print $6}')
    liftover_to_hg38=$(echo "$row" | awk -F'\t' '{print $7}')
    liftover_to_mm10=$(echo "$row" | awk -F'\t' '{print $8}')
    
    # Combine non-empty chain URLs with commas
    chain_urls=""
    for url in "$liftover_from_hg38" "$liftover_from_mm10" "$liftover_to_hg38" "$liftover_to_mm10"; do
        if [[ -n "$url" ]]; then
            if [[ -n "$chain_urls" ]]; then
                chain_urls="${chain_urls},${url}"
            else
                chain_urls="$url"
            fi
        fi
    done
    
    # Submit job with chain URLs if any exist
    echo "Submitting job for target genome: $genome"
    if [[ -n "$chain_urls" ]]; then
        sbatch --job-name="dl_${genome}" scripts/download-genome.sh \
            -g "$genome" \
            -f "$fasta_url" \
            -p "$PROJECT_DIR" \
            -s "$SCRATCH_DIR" \
            -l "$chain_urls"
    else
        sbatch --job-name="dl_${genome}" scripts/download-genome.sh \
            -g "$genome" \
            -f "$fasta_url" \
            -p "$PROJECT_DIR" \
            -s "$SCRATCH_DIR"
    fi
done

echo "All download jobs submitted"
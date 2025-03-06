#!/bin/bash
#
# Check and generate liftover chains for target genomes
#
PROJECT_DIR=$HOME/repos/custom_ArchR_genomes_and_annotations
SCRATCH_DIR=/scratch/$USER
GENOME_DIR=${PROJECT_DIR}/genomes

# Process genomes from config file
conda activate custom_genes
echo "Checking for missing liftover chains..."

config_file="config/target_genomes.tsv"
num_genomes=$(wc -l < $config_file)

# Process each genome starting from line 2 (skip header)
for i in $(seq 2 $num_genomes); do
    # Parse TSV line using awk
    row=$(awk -v line=$i 'NR==line' "$config_file")
    genome=$(echo "$row" | awk -F'\t' '{print $2}')
    chain_hg38=$(echo "$row" | awk -F'\t' '{print $5}')
    chain_mm10=$(echo "$row" | awk -F'\t' '{print $6}')
    
    # Skip if chains already exist
    if [[ -n "$chain_hg38" || -n "$chain_mm10" ]]; then
        echo "Chains already exist for $genome, skipping..."
        continue
    fi
    
    echo "No chains found for $genome, submitting jobs..."
    
    # Set up paths
    source_fasta="${GENOME_DIR}/${genome}/${genome}.fa.gz"
    hg38_fasta="${GENOME_DIR}/hg38/hg38.fa.gz"
    mm10_fasta="${GENOME_DIR}/mm10/mm10.fa.gz"
    chain_dir="${GENOME_DIR}/${genome}/liftover"
    mkdir -p "$chain_dir"
    
    # Determine proximity based on genome name pattern
    if [[ $genome =~ ^(rheMac|macFas|mCalJac|calJac|mMacNem) ]]; then
        proximity="primates"
    else
        proximity="mammals"
    fi
    
    # Submit chain generation jobs
    sbatch --job-name="chain_${genome}_hg38" scripts/make-liftover-chain.sh \
        -t "$source_fasta" \
        -q "$hg38_fasta" \
        -p "$proximity" \
        -o "$chain_dir" \
        -w "$SCRATCH_DIR"
           
    sbatch --job-name="chain_${genome}_mm10" scripts/make-liftover-chain.sh \
        -t "$source_fasta" \
        -q "$mm10_fasta" \
        -p "$proximity" \
        -o "$chain_dir" \
        -w "$SCRATCH_DIR"
done

echo "All liftover chain jobs submitted"



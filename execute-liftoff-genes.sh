#!/bin/bash
#
# Execute liftoff jobs for all source/target genome combinations
#

PROJECT_DIR=${1:-$HOME/repos/custom_ArchR_genomes_and_annotations}
SCRATCH_DIR=${2:-/scratch/$USER}

echo "Project directory: $PROJECT_DIR"
echo "Scratch directory: $SCRATCH_DIR"

# Process source genomes config
echo "Reading source genomes config..."
config_file="$PROJECT_DIR/config/source_genomes.tsv"
num_sources=$(wc -l < "$config_file")

# Get list of target genomes
echo "Reading target genomes config..."
targets_file="$PROJECT_DIR/config/target_genomes.tsv"
num_targets=$(wc -l < "$targets_file")

# Create logs directory
mkdir -p "$PROJECT_DIR/logs"

# Loop through human genomes to non-rodents
for i in $(seq 2 5); do
    # Parse source genome info
    row=$(awk -v line=$i 'NR==line' "$config_file")
    source_genome=$(echo "$row" | awk -F'\t' '{print $2}')
    annot_file=$(echo "$row" | awk -F'\t' '{print $5}' | tr -d '\r\n')
    annot_name=$(echo "$row" | awk -F'\t' '{print $5}' | tr -d '\r\n')
    
    # Loop through target genomes
    # for j in $(seq 2 $num_targets); do
    for j in $(seq 2 8); do
        # Parse target genome info
        row=$(awk -v line=$j 'NR==line' "$targets_file")
        target_genome=$(echo "$row" | awk -F'\t' '{print $2}' | tr -d '\r\n')
        
        echo "Submitting liftoff job for: $source_genome -> $target_genome"
        
        # Submit SLURM job
        sbatch --job-name="liftoff_${source_genome}_${target_genome}" \
            "$PROJECT_DIR/scripts/liftoff-genes.sh" \
            -s "$source_genome" \
            -t "$target_genome" \
            -a "$annot_name" \
            -p "$PROJECT_DIR" \
            -w "$SCRATCH_DIR"
    done
done

# for the mouse gene annotations to rats
for i in $(seq 6 8); do
    # Parse source genome info
    row=$(awk -v line=$i 'NR==line' "$config_file")
    source_genome=$(echo "$row" | awk -F'\t' '{print $2}')
    annot_file=$(echo "$row" | awk -F'\t' '{print $5}' | tr -d '\r\n')
    annot_name=$(echo "$row" | awk -F'\t' '{print $5}' | tr -d '\r\n')
    
    # Loop through target genomes
    # for j in $(seq 2 $num_targets); do
    for j in $(seq 9 10); do
        # Parse target genome info
        row=$(awk -v line=$j 'NR==line' "$targets_file")
        target_genome=$(echo "$row" | awk -F'\t' '{print $2}' | tr -d '\r\n')
        
        echo "Submitting liftoff job for: $source_genome -> $target_genome"
        
        # Submit SLURM job
        sbatch --job-name="liftoff_${source_genome}_${target_genome}" \
            "$PROJECT_DIR/scripts/liftoff-genes.sh" \
            -s "$source_genome" \
            -t "$target_genome" \
            -a "$annot_name" \
            -p "$PROJECT_DIR" \
            -w "$SCRATCH_DIR"
    done
done

echo "All liftoff jobs submitted"
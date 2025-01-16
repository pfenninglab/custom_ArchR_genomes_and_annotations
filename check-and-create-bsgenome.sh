#!/bin/bash
#
# Execute BSgenome builds from config file
#


PROJECT_DIR=$HOME/repos/custom_ArchR_genomes_and_annotations
SCRATCH_DIR=/scratch/$USER
GENOME_DIR=${PROJECT_DIR}/output/genomes

# Process genomes from config file
conda activate custom_genes

echo "Processing genomes..."
config_file="config/target_genomes.tsv"
num_genomes=$(wc -l < $config_file)

# Process each genome starting from line 2 (skip header)
for i in $(seq 2 $num_genomes); do
    # Parse TSV line using awk
    row=$(awk -v line=$i 'NR==line' "$config_file")
    genome=$(echo "$row" | awk -F'\t' '{print $2}')
    species=$(echo "$row" | awk -F'\t' '{print $3}')
    source_file=$(echo "$row" | awk -F'\t' '{print $4}')
    pgk_name=$(echo "$row" | awk -F'\t' '{print $5}')
    fasta_path=${GENOME_DIR}/${genome}/${genome}.fa.gz
    
    if [ -z "$pgk_name" ]; then
        # Time to create the package name
        GENUS=$(echo "$species" | cut -d' ' -f1)
        SPECIES_NAME=$(echo "$species" | cut -d' ' -f2)
        SPECIES_FORMATTED=$(echo "${GENUS:0:1}${SPECIES_NAME}" | sed 's/^./\U&/g; s/ //g')
        pgk_name="BSgenome.${SPECIES_FORMATTED}.Custom.${genome}"
    else
        continue
    fi
    
    if [[ "$pkg_name" == *"Custom"* ]]; then
        if Rscript -e "if('${pgk_name}' %in% rownames(installed.packages())) quit(status = 0) else quit(status = 1)"; then
            echo "Package ${pgk_name} already installed, skipping..."
        else
            echo "Submitting job for genome: $genome"
            # Uncomment the line below to actually submit the job
            sbatch --job-name="bsg_${genome}" scripts/make-bsgenome.sh \
                -s "$species" \
                -g "$genome" \
                -f "$fasta_path" \
                -n "$pgk_name" \
                -w "$SCRATCH_DIR" \
                -o "${GENOME_DIR}/${genome}" \
                -p "$PROJECT_DIR"
        fi
    fi
    
done

echo "All BSgenome build jobs submitted"
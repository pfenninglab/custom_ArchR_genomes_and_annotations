# this script is to borrow the mitochondrial chromosomes from earlier genome versions of certain species to append to newer genome versions

PROJDIR=/home/bnphan/repos/custom_ArchR_genomes_and_annotations
DATADIR=${PROJDIR}/genomes
SCRATCH=/scratch/${USER}

mkdir -p $SCRATCH

borrow_mitochondrial_chromosome() {
    # Parse arguments
    while getopts s:t:m:d: flag; do
        case "${flag}" in
            s) SOURCE_GENOME_GZ=${OPTARG};;
            t) TARGET_GENOME_GZ=${OPTARG};;
            m) MT_CHR_NAME=${OPTARG};;
            d) OUTPUT_DIR=${OPTARG};;
        esac
    done

    # Validate required parameters
    if [[ -z "$SOURCE_GENOME_GZ" || -z "$TARGET_GENOME_GZ" || -z "$MT_CHR_NAME" || -z "$OUTPUT_DIR" ]]; then
        echo "Error: Missing required parameters"
        echo "Usage: borrow_mitochondrial_chromosome -s <source_genome.gz> -t <target_genome.gz> -m <mt_chr_name> -d <output_dir>"
        return 1
    fi

    # Check if input files exist
    if [[ ! -f "$SOURCE_GENOME_GZ" ]]; then
        echo "Error: Source genome file not found: $SOURCE_GENOME_GZ"
        return 1
    fi
    if [[ ! -f "$TARGET_GENOME_GZ" ]]; then
        echo "Error: Target genome file not found: $TARGET_GENOME_GZ"
        return 1
    fi

    # Set up paths and filenames
    SCRATCH="/scratch/${USER}"
    mkdir -p "$SCRATCH"
    mkdir -p "$OUTPUT_DIR"

    SOURCE_BASENAME=$(basename "$SOURCE_GENOME_GZ" .fa.gz)
    TARGET_BASENAME=$(basename "$TARGET_GENOME_GZ" .fa.gz)
    FINAL_OUTPUT="${OUTPUT_DIR}/${TARGET_BASENAME}.plusMT.fa.gz"

    # Check if final output already exists
    if [[ -f "$FINAL_OUTPUT" ]]; then
        echo "Output file already exists: $FINAL_OUTPUT"
        echo "Use a different output directory or remove the existing file to proceed."
        return 0
    fi

    # Prepare scratch file paths
    SCRATCH_SOURCE_GENOME_GZ="${SCRATCH}/${SOURCE_BASENAME}.fa.gz"
    SCRATCH_TARGET_GENOME_GZ="${SCRATCH}/${TARGET_BASENAME}.fa.gz"
    SOURCE_GENOME="${SCRATCH}/${SOURCE_BASENAME}.fa"
    TARGET_GENOME="${SCRATCH}/${TARGET_BASENAME}.fa"
    MT_CHROMOSOME="${SCRATCH}/${SOURCE_BASENAME}.${MT_CHR_NAME}.fa"
    TARGET_GENOME_MT="${SCRATCH}/${TARGET_BASENAME}.plusMT.fa"

    # Copy gzipped genomes to scratch with progress indicators
    echo "Copying source genome to scratch..."
    rsync -Pa "$SOURCE_GENOME_GZ" "$SCRATCH_SOURCE_GENOME_GZ" || { echo "Error: Failed to copy source genome"; return 1; }
    
    echo "Copying target genome to scratch..."
    rsync -Pa "$TARGET_GENOME_GZ" "$SCRATCH_TARGET_GENOME_GZ" || { echo "Error: Failed to copy target genome"; return 1; }

    # Extract genomes only if needed
    if [[ ! -f "$SOURCE_GENOME" ]]; then
        echo "Extracting source genome..."
        gunzip -c "$SCRATCH_SOURCE_GENOME_GZ" > "$SOURCE_GENOME" || { echo "Error: Failed to extract source genome"; return 1; }
    fi
    
    if [[ ! -f "$TARGET_GENOME" ]]; then
        echo "Extracting target genome..."
        gunzip -c "$SCRATCH_TARGET_GENOME_GZ" > "$TARGET_GENOME" || { echo "Error: Failed to extract target genome"; return 1; }
    fi

    # Extract mitochondrial chromosome from the source species
    echo "Extracting mitochondrial chromosome: $MT_CHR_NAME"
    if [[ ! -f "${SOURCE_GENOME}.fai" ]]; then
        samtools faidx "$SOURCE_GENOME" || { echo "Error: Failed to index source genome"; return 1; }
    fi
    
    # Check if the MT chromosome exists in the source genome
    if ! grep -q "^$MT_CHR_NAME" "${SOURCE_GENOME}.fai"; then
        echo "Error: Mitochondrial chromosome '$MT_CHR_NAME' not found in source genome"
        return 1
    fi
    
    samtools faidx "$SOURCE_GENOME" "$MT_CHR_NAME" > "$MT_CHROMOSOME" || { 
        echo "Error: Failed to extract mitochondrial chromosome"; 
        return 1; 
    }

    # Append mitochondrial chromosome to target genome
    echo "Creating combined genome..."
    cat "$TARGET_GENOME" "$MT_CHROMOSOME" > "$TARGET_GENOME_MT" || { 
        echo "Error: Failed to create combined genome"; 
        return 1; 
    }
    
    echo "Indexing combined genome..."
    samtools faidx "$TARGET_GENOME_MT" || { echo "Error: Failed to index combined genome"; return 1; }

    # Compress final genome and move to output directory
    echo "Compressing final genome..."
    gzip -c "$TARGET_GENOME_MT" > "${TARGET_GENOME_MT}.gz" || { 
        echo "Error: Failed to compress combined genome"; 
        return 1; 
    }
    
    echo "Moving final genome to output directory..."
    rsync -Pa "${TARGET_GENOME_MT}.gz" "$FINAL_OUTPUT" || { 
        echo "Error: Failed to copy final genome to output directory"; 
        return 1; 
    }

    # Clean up temporary files (optional, uncomment if desired)
    # echo "Cleaning up temporary files..."
    # rm -f "$SOURCE_GENOME" "$TARGET_GENOME" "$MT_CHROMOSOME" "$TARGET_GENOME_MT" "${TARGET_GENOME_MT}.gz" "${TARGET_GENOME_MT}.fai"

    echo "Success: Mitochondrial chromosome (${MT_CHR_NAME}) added to $FINAL_OUTPUT"
    return 0
}
# Example usage:
# borrow_mitochondrial_chromosome -s /path/to/macFas5.fa.gz -t /path/to/macFas6.fa.gz -m chrM -d /path

###################################
# for the Macaca Fasicularis
MT_CHR_NAME='chrM'
SOURCE_SPECIES=macFas5
TARGET_SPECIES=macFas6

borrow_mitochondrial_chromosome \
-s ${DATADIR}/${SOURCE_SPECIES}/${SOURCE_SPECIES}.fa.gz \
-t ${DATADIR}/${TARGET_SPECIES}/${TARGET_SPECIES}.fa.gz \
-m $MT_CHR_NAME -d ${DATADIR}/${TARGET_SPECIES}

mv ${DATADIR}/${TARGET_SPECIES} ${DATADIR}/${TARGET_SPECIES}.plusMT

###################################
# for the marmosets
MT_CHR_NAME='chrM'
SOURCE_SPECIES=calJac4
TARGET_SPECIES=mCalJa1.2.pat.X

borrow_mitochondrial_chromosome \
-s ${DATADIR}/${SOURCE_SPECIES}/${SOURCE_SPECIES}.fa.gz \
-t ${DATADIR}/${TARGET_SPECIES}/${TARGET_SPECIES}.fa.gz \
-m $MT_CHR_NAME -d ${DATADIR}/${TARGET_SPECIES}

mv ${DATADIR}/${TARGET_SPECIES} ${DATADIR}/${TARGET_SPECIES}.plusMT




# Gene Annotation Lift - Automated Cross-Species Gene Annotation Pipeline

## Overview

This pipeline automates the process of mapping high-quality gene annotations from well-annotated reference genomes (human/mouse) to other species. It generates:

1. Mapped gene annotations for target species using liftoff
2. ArchR-compatible genome and gene annotation files
3. Quality control metrics and comparisons between different annotation sources

## Key Features

- Flexible source genome selection (human/mouse versions)
- Configurable target species and genome versions
- Automated generation of ArchR-compatible files
- Progress tracking and intelligent re-running
- Quality control and annotation comparison metrics

## Directory Structure

```
.
├── config/
│   ├── source_genomes.tsv      # Reference genome configs 
│   ├── target_genomes.tsv      # Target species configs
│   └── run_config.yaml         # Pipeline parameters
├── scripts/
│   ├── 01_download_genomes.sh  # Download source/target genomes
│   ├── 02_liftoff_genes.sh     # Run liftoff mapping
│   ├── 03_make_archr.R         # Generate ArchR files
│   └── 04_qc_annotations.R     # Compare annotation quality
├── output/
│   ├── genomes/               # Downloaded genome files
│   ├── annotations/           # Mapped gene annotations
│   ├── archr/                 # ArchR-compatible files  
│   └── qc/                    # Quality metrics
└── logs/                      # Pipeline run logs
```

## Configuration

### Source Genomes (source_genomes.tsv)
```tsv
source_id  species  genome_ver  annotation_ver  annotation_url
hg38       human    GRCh38     GENCODE_v38     ftp://...
mm10       mouse    GRCm38     GENCODE_vM25    ftp://...
```

### Target Genomes (target_genomes.tsv) 
```tsv
target_id     species            genome_ver  genome_url
rheMac10      rhesus_macaque    rheMac10    ftp://...
mCalJac1      marmoset          mCalJac1    https://...
rn7           rat               mRatBN7.2   ftp://...
```

## Usage

1. Edit configuration files in `config/`

2. Run the full pipeline:
```bash
snakemake --configfile config/run_config.yaml
```

Or run individual steps:
```bash
# Download genomes
bash scripts/01_download_genomes.sh

# Run liftoff gene mapping 
bash scripts/02_liftoff_genes.sh

# Generate ArchR files
Rscript scripts/03_make_archr.R

# Compare annotations
Rscript scripts/04_qc_annotations.R
```

## Output Files

For each target species, the pipeline generates:

1. Gene Annotations:
- `{target_id}_liftoff_{source_id}.gff3` - Mapped annotations in GFF3 format
- `{target_id}_liftoff_{source_id}.gtf` - Mapped annotations in GTF format

2. ArchR Files:
- `{target_id}_ArchR_annotations.rda` - ArchR genome and gene annotations
- `{target_id}_genes_exons_TSS.rda` - Gene features for ArchR

3. QC Metrics:
- `{target_id}_annotation_comparisons.pdf` - Comparison visualizations
- `{target_id}_mapping_stats.tsv` - Mapping quality metrics

## Creating an R Package

After running the pipeline, use the following script to package the outputs into a distributable R package:

```R
source("scripts/create_package.R")
```

This will:
1. Create proper R package structure
2. Add documentation for genome and gene annotations
3. Package ArchR-compatible files
4. Generate installation instructions

## Citation

If you use this pipeline in your research, please cite:

[Citation information]

## Contributing

We welcome contributions! Please see CONTRIBUTING.md for guidelines.

## License

[License information]
# Custom genome and gene annotations for cross-species genomics

This repository provides a systematic framework for mapping gene annotations between species, with a focus on mapping human (GENCODE v47) and mouse (GENCODE vM25) annotations to various target species including primates and rodents.

## Overview

Cross-species genomics requires high-quality genome assemblies and gene annotations. While genome assemblies are increasingly available through efforts like the [Vertebrate Genome Project](https://vertebrategenomesproject.org), gene annotations often lag behind. This repository provides:

1. Automated downloading of source and target genome assemblies
2. Systematic mapping of GENCODE annotations using [Liftoff](https://github.com/agshumate/Liftoff)
3. Creation of genome packages for downstream analysis

## Source Data

### Source Genomes and Annotations
Links to genome sequence and gene annotations for human and mouse also located in `config/source_genomes.tsv`:


| Genome | Version | Source | Annotation |
|--------|---------|--------|------------|
| Human | GRCh38/hg38 | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) | [GENCODE v47 comprehensive](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz) |
| Human | GRCh38/hg38 | UCSC | [GENCODE v47 basic](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz) |
| Mouse | GRCm38/mm10 | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz) | [GENCODE vM25](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz) |

### Target Genomes  
Links to genome sequence for target genomes also located in `config/target_genomes.tsv`:


- Rhesus macaque
  - rheMac8 ([UCSC download](https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/rheMac8.fa.gz))
  - rheMac10 ([UCSC download](https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz))
- Crab-eating macaque
  - macFas6 ([NCBI download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/100/615/GCA_011100615.1_Macaca_fascicularis_6.0/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.gz))
- Marmoset
  - mCalJac1 ([VGP Genomeark download](https://s3.amazonaws.com/genomeark/species/Callithrix_jacchus/mCalJac1/assembly_curated/mCalJac1.mat.cur.20200212.fasta.gz))
  - calJac4 ([UCSC download](https://hgdownload.soe.ucsc.edu/goldenPath/calJac4/bigZips/calJac4.fa.gz))
- Pig-tailed macaque
  - mMacNem1 ([VGP Genomeark download](https://s3.amazonaws.com/genomeark/species/Macaca_nemestrina/mMacNem1/assembly_curated/mMacNem1.hap1.cur.20240610.fasta.gz))
- Norway rat
  - rn6 ([UCSC download](https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz))
  - rn7 ([UCSC download](https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.gz))
- Pig
  - susScr11 ([UCSC download](https://hgdownload.soe.ucsc.edu/goldenPath/susScr11/bigZips/susScr11.fa.gz))

## Output Structure

Processed data is organized under `output/genomes/` with the following structure:
The download utility below uses data from the source and target genome .tsv files
to download the required genome and gene annotations for the runs.
```
output/genomes/
├── {genome}/              # e.g., rheMac10/
│   ├── {genome}.fa       # Genome FASTA
│   └── annotations/      
│       └── {genome}-{source}-{version}.gtf.gz  # Lifted annotations
```

Running the liftoff wrapper will also create the mapping from human/mouse to the 
target genomes. 

Example for rheMac10:
```
rheMac10/
├── rheMac10.fa
└── annotations/
    ├── rheMac10-hg38-gencode.v47.basic.gtf.gz
    ├── rheMac10-hg38-gencode.v47.comp.gtf.gz
    └── rheMac10-mm10-gencode.vM25.gtf.gz
```

## Usage

### Requirements
- Install dependencies using provided conda environment:
```bash
# clone the repository
git clone git@github.com:pfenninglab/custom_ArchR_genomes_and_annotations.git
cd custom_ArchR_genomes_and_annotations

# create conda environment from config file
conda env create -f config/conda_environment.yml
```

### Basic Workflow

1. Download source/target genomes:
```bash
bash scripts/download-genome.sh -g rheMac10 \
  -f https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz \
  -n gencode.v47.basic \
  -t https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz
```

2. Run Liftoff gene mapping:
```bash
bash scripts/liftoff-genes.sh -s hg38 -t rheMac10 -a gencode.v47.basic
```

## Citations

If you use these resources, please cite:

- Liftoff: [Shumate A, et al. (2021) Bioinformatics](https://academic.oup.com/bioinformatics/article/37/12/1639/6035128)
- GENCODE: [Frankish A, et al. (2021) Nucleic Acids Research](https://academic.oup.com/nar/article/49/D1/D916/6018430)
- Original gene annotations: Please cite the data resource DOI:
```
Phan, BaDoi; Pfenning, Andreas (2022): Alternate gene annotations for rat, macaque, and marmoset for single cell RNA and ATAC analyses.
Carnegie Mellon University. Dataset. https://doi.org/10.1184/R1/21176401.v1
```

## Internal Usage
If you are a Pfenning lab member with access to the Lane cluster, this repository
and all final github files and intermediate files are located at
```bash
cd /home/bnphan/repos/custom_ArchR_genomes_and_annotations
```

## Contributing

Issues and pull requests welcome! See CONTRIBUTING.md for guidelines.
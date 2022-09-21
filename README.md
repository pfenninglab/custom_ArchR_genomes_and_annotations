# Custom ArchR genome and gene annotations
## by BaDoi Phan (badoi dot phan at pitt dot edu)

# 1) Premise: 
Not all of single-cell ATAC-seq biomedical molecular epigenetics is done in human and mouse genomes where there are 
high quality genomes and gene annotations. For the other species that are still highly relevant to study health and disease, 
here are some ArchR annotations to enable less frustration to have snATAC-seq data analyzed with [ArchR](https://www.archrproject.com). 

# 2) Strategy for gene annotations:
We can use the evolution of related mammalian species tend to have orthologous gene elements (TSS, exons, genes). For example, house mouse (**mus musculus**) is a median of 15.4MY diverged from the Norway rat (_rattus norvegicus_), with [TimeTree](http://www.timetree.org). Humans are a median of 28.9 MY diverged from rhesus macaques. To borrow the higher quality and more complete gene annotations, we can use a gene-aware method of lifting gene annotations from one genome to another, [liftoff, Shumate and Salzberg, 2021](https://academic.oup.com/bioinformatics/article/37/12/1639/6035128). For the source of "high quality" gene annotation, we use the NCBI Refseq annotations from the *hg38/GRCh38* and *mm10/GRCm38* annotations downloaded from the UCSC Genome browser. 

Similarly, much work by the ENCODE Consortium has found with the large human and mouse epigenomic data that certain regions of the genome in these species have artifactual signals and need to be excluded from epigenomic analsyes, [Amemiya et al., 2021](https://www.nature.com/articles/s41598-019-45839-z). These regions were pulled from and human and mouse from [here](https://github.com/Boyle-Lab/Blacklist/) and used the liftOver to map to the target genomes below, for simplicity.

# 3) list of resources by file name
Surprisingly, all these files are small enough to put on github for a couple custom genomes. Below are the organizations 
- \*liftoff\*.gtf.gz and \*liftoff\*.gff3.gz: the gzipped annotation from the higher quality annotations to the target genome using [liftoff](https://github.com/agshumate/Liftoff)
- \*liftOver\*blacklist.v2.bed: the ENCODE regions to exclude from epigenomic analyses mapped to the target genome using [liftOver](https://genome-store.ucsc.edu)
- \*ArchRGenome.R: the Rscript used to make the custom ArchR annotations
- \*ArchR_annotations.rda: the R Data object that contains the *geneAnnotation* and *genomeAnnotation* objects to use with [ArchR::createArrowFiles()](https://www.archrproject.com/reference/createArrowFiles.html)

# 4) list of species/genomes/source files
For most of these files, the genome fasta sequences were grabbed from the UCSC Genome Browser at https://hgdownload.soe.ucsc.edu/goldenPath/${GENOME_VERSION}/, where ${GENOME_VERSION} is any of the version below except **mCalJac1**. Some of these genomes were updated from the Vertebrate Genome Project, which seeks to create complete rather than draft genome assemblies of all mammals on the planet, [Rhie et al. 2021](https://www.nature.com/articles/s41586-021-03451-0). These genomes have **VGP** and that naming version if there's an alternate naming scheme. The VGP is [pretty cool](https://vertebrategenomesproject.org) and they make good [genome assemblies](https://vgp.github.io/genomeark/).

- rn6: [rat genome v6, BCM-Baylor version](https://www.nature.com/articles/nature02426)
- rn7: [rat genome also called VGP mRatBN7.2](https://journals.physiology.org/doi/abs/10.1152/physiolgenomics.00017.2022)
- rheMac8: [rhesus macaque v8](https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/)
- rheMac10: [rhesus macaque v10](https://www.science.org/doi/10.1126/science.abc6617?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
- mCalJac1: marmoset VGP genome, [fasta from the maternal assembly here](https://www.ncbi.nlm.nih.gov/assembly/GCA_011078405.1/)


## Internal
These files are a subset of resources at the CMU Computational Biology Department Lane Cluster at:
`/home/bnphan/resources/genomes`

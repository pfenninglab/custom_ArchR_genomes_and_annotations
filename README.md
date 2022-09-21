# Custom ArchR genome and gene annotations
## by BaDoi Phan (badoi dot phan at pitt dot edu)

# 1) Premise: 
Not all of single-cell ATAC-seq biomedical molecular epigenetics is done in human and mouse genomes where there are 
high quality genomes and gene annotations. For the other species that are still highly relevant to study health and disease, 
here are some ArchR annotations to enable less frustration to have snATAC-seq data analyzed with [ArchR](https://www.archrproject.com). 

# 2) Strategy for gene annotations:
We can use the evolution of related mammalian species tend to have orthologous gene elements (TSS, exons, genes). For example, house mouse (**mus musculus**) is a median of 15.4MY diverged from the Norway rat (_rattus norvegicus_), with (TimeTree)[http://www.timetree.org]. Humans are a median of 28.9 MY diverged from rhesus macaques. To borrow the higher quality and more complete gene annotations, we can use a gene-aware method of lifting gene annotations from one genome to another, (liftoff, Shumate and Salzberg, 2021)[https://academic.oup.com/bioinformatics/article/37/12/1639/6035128]. For the source of "high quality" gene annotation, we use the NCBI Refseq annotations from the *hg38/GRCh38* and *mm10/GRCm38* annotations downloaded from the UCSC Genome browser.

# 3) list of resources by file name
- *.R: likely the Rscript used to make the custom ArchR *geneAnnotation* and *genomeAnnotation* objects to use with (ArchR::createArrowFiles())[https://www.archrproject.com/reference/createArrowFiles.html]


# 4) list of species/genomes



## Default repo
local({r <- getOption("repos")
    r["CRAN"] <- "http://cran.r-project.org" 
    options(repos=r)
})

#this script is meant to be run in a clean R environment, after having created a conda environment from ArchR_Env.yaml
#packages to install with R's default install command

normalPackages = c("BiocManager", "pheatmap", "hexbin", "rhandsontable")
install.packages(setdiff(normalPackages, rownames(installed.packages()))) 

#packages to be installed with devtools
library("stringr")
devPackages = c("immunogenomics/presto", "immunogenomics/harmony", "GreenleafLab/chromVARmotifs")
devPackagesNames = gsub("/", "", str_extract(devPackages, "/.*"))
if (!require(devPackagesNames, quietly = TRUE))
    devtools::install_github(devPackages)  

#Seurat depends on normal and devtools packages
install.packages(setdiff("Seurat", rownames(installed.packages())))

#packages to be installed with BioCManager
bioCpackages = c("rtracklayer", "CNEr", "TFBSTools", "motifmatchr", "restfulr", "annotate", "BSGenome", "chromVAR")
BiocManager::install(setdiff(bioCpackages, rownames(installed.packages())))


#install ArchR
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library()

R.version

installed = dir(.libPaths()[1])
installed

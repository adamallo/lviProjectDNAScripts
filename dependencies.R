#!/usr/bin/env Rscript

#CRAN
requiredCranPackages <- c("devtools",
                      "BiocManager",
                      "tidyverse",
                      "ape",
                      "cowplot",
                      "data.table",
                      "emmeans",
                      "ggbeeswarm",
                      "ggpubr",
                      "ggsignif",
                      "ggtree",
                      "ggimage",
                      "glmmTMB",
                      "matrixStats",
                      "parallel",
                      "phangorn",
                      "tidytree")

toInstallCranPackages <- requiredCranPackages[!requiredCranPackages %in% installed.packages()]
install.packages(toInstallCranPackages, dependencies = T, repos = "http://cran.us.r-project.org")

#Bioconductor
library("BiocManager")
requiredBiocPackages <- c("aCGH",
                          "AnnotationDbi",
                          "CGHcall",
                          "HDO.db",
                          "clusterProfiler",
                          "DNAcopy",
                          "enrichplot",
                          "GenomicFeatures",
                          "GenomicRanges",
                          "org.Hs.eg.db",
                          "QDNAseq",
                          "rtracklayer",
                          "TxDb.Hsapiens.UCSC.hg38.knownGene")

toInstallBiocPackages <- requiredBiocPackages[!requiredBiocPackages %in% installed.packages()]
BiocManager::install(toInstallBiocPackages,update = T,ask = F)

#Github packages
library("devtools")
requiredGithubPackages <- c("stephenturner/annotables",
                            "Sawyer-s-Group/breakclone",
                            "crukci-bioinformatics/rascal")
toInstallGithubPackages <- requiredGithubPackages[!gsub("[^/]*/","",requiredGithubPackages) %in% installed.packages()]
lapply(toInstallGithubPackages,function(package){install_github(package,dependencies = T,upgrade = T)})


# README
Repository to replicate the DNA-seq and downstream analyses in: Characterization of the lymphovascular invasion microenvironment reveals immune response dichotomy. Rivero-Gutierrez, Mallo, et al. 2024

# Description
These analysis are carried out in eight main steps, CNA calling, Genetic analyses, Clonality analyses, recurrent CNA calling, analysis of the associations between CNAs and stages and subtypes, comparison between differentially-expressed and differentially-altered genes, copy number alteration enrichment analyses, and phylogenetic reconstruction. CNA calling is a pre-requisite for all of them.

# Instructions
To replicate these analyses, follow these steps:
1. Modify the *configFile.in* file to configure it for your system and rename it as *configFile*. You can find an example configuration file in the repository, *configFile.ex*.
2. Export an environmental variable with the absolute path of this repository `export lviProjectDNAScripts="/path/to/the/repository/lviProjectDNAScripts/"`
3. Make sure to have installed all the dependencies [listed below](#dependencies). You can execute the R script *dependencies.R* to install all R packages.
4. Run the analyses you are interested in from the list below

## Analyses
- **CNA calling**: the script *CNACalling.R* performs the CNA calling using QDNAseq and Rascal. This step is required for all downstream analyses.

- **Genetic analyses**: the script *genetics.R* performs the statistical analysis and plots for the genetic characterization of the different neoplasia stages.

- **Clonality analysis**: the script *clonalityScore.R* assesses the clonality of the different neoplastic stages using breakclone.

- **GISTIC 2.0 analysis**: the script *makeGisticInputFromTsvFiles.sh* puts together the input file for GISTIC 2.0.
    1. Execute *makeGisticInputFromTsvFiles.sh* as `source $lviProjectDNAScripts/configFile; $scriptsDir/makeGisticInputFromTsvFiles.sh -o $gisticInputFile $(find $acnaDir -name "*highcell*.tsv" | grep -v NORMAL)`
    2. Run GISTIC2.0 using the following command `source $lviProjectDNAScripts/configFile; gistic2 -b $gisticOutDir -seg $gisticInputFile -refgene $gisticMatrixFile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme`.
    - Alternatively, we also used only one LVI sample per patient instead of all non-normal samples to detect recurrent alterations in LVI tissue only. For this, run this alternative series of steps after step 1.:
    2. Execute *lviSelection.R*
    3. Run GISTIC2.0 using the following command `source $lviProjectDNAScripts/configFile; gistic2 -b $gisticLVIOutDir -seg $acnaDir/$gisticLVIInputFile -refgene $gisticMatrixFile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme`.

- **Associations between recurrent CNAs and Stages and LVI subtypes**: the script *gisticAssociations.R* does these comparisons and requires the GISTIC 2.0 results.

- **Comparison between differential expression and copy number alteration results**: the script *diffExpVsCNAs.R* does these comparisons.

- **Enrichment analyses**: CNA enrichment analyses of all alterations (*deNovoEnrichment.R*) or recurrent alterations (*deNovoEnrichmentGistic.R*). The second requires the GISTIC 2.results.

- **Phylogenetic/dendogram reconstruction**: the scripts *getBinaryData.R* generates binary data using breakpoint information, which then is used in Phylip to reconstruct the best sample tree and its bootstrap support (*runPhylip.sh*, must be run from the directory with the results from *getBinaryData.R*, for example using `$scriptsDir/runPhylip.sh`. It uses *dollopCMD*, *dollopCMDBoost*, and *seqbootCMD* command files internally). Finally, the script *dollopAnalysis.R* performs the ancestral state reconstruction and generates the plots for the manuscript.

## Dependencies

### Other software
- Perl
- GISTIC2.0
- Phylip (you can find a version with minor makefile modifications for compilation on Apple Silicon systems and pre-compiled binary versions at https://github.com/adamallo/phylip)

### R packages
- aCGH
- annotables #devtools::install_github("stephenturner/annotables")
- AnnotationDbi
- ape
- breakclone ##devtools::install_github("Sawyer-s-Group/breakclone")
- car
- CGHcall
- clusterProfiler
- cowplot
- devtools
- data.table
- DNAcopy
- dplyr
- emmeans
- enrichplot
- GenomicFeatures
- GenomicRanges
- ggbeeswarm
- ggplot2
- ggpubr
- ggsignif
- ggtree
- glmmTMB
- matrixStats
- org.Hs.eg.db
- parallel
- phangorn
- QDNAseq
- rascal ##devtools::install_github("crukci-bioinformatics/rascal")
- readxl
- rtracklayer
- tibble
- tidytree
- tidyverse
- TxDb.Hsapiens.UCSC.hg38.knownGene

## Data:
Data needed to run these analyses:
- rawReads: Read counts for 50Kbp bins
- designFile.csv: File with sample information
- diffExpFiles: Files with differentially-expressed genes

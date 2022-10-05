README

Repository to replicate the DNA-seq and downstream analyses in: Characterization of the lymphovascular invasion microenvironment reveals immune response dichotomy. Rivero-Gutierrez, Mallo, et al. 2022

These analysis are carried out in three main steps:

- **CNA calling**: the script CNACalling.R performs the CNA calling using QDNAseq and Rascal. This step is required for both the phylogenetic and genetic analyses.

- **Genetic analyses**: the script genetics.R performs the statistical analysis and plots for the genetic characterization of the different neoplasia stages. Requires the output of CNACalling.R

- **Phylogenetic/dendogram reconstruction**: the scripts getBinaryData.R generates binary data using breakpoint information, which then is used in Phylip to reconstruct the best sample tree and its bootstrap support (commands in runPhylip.sh, run manually). Finally, the script dollopAnalysis.R performs the ancestral state reconstruction and generates the plots for the manuscript. dollopCMD, dollopCMDBoost, and seqbootCMD 

Data (pending to be added to this repository or deposited somewhere else):

- rawReads: Read counts for 50Kbp bins
- designFile.csv: File with sample information

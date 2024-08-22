suppressMessages({
  library(tidyverse)
  library(readxl)
  library(parallel)
  library(annotables) #devtools::install_github("stephenturner/annotables")
})

#FUNCTIONS
calcOverlap <- function(refStart,refEnd,segStart,segEnd){
  overlapStart=ifelse(refStart>segStart,refStart,segStart)
  overlapEnd=ifelse(refEnd<segEnd,refEnd,segEnd)
  return(overlapEnd-overlapStart)
}

wilcox.test.wErrorHandling <- function(x,y,...){
  tryCatch({
    withCallingHandlers(
      {
        wilcox.test(x,y,...)
      },
      warning = function(w){
        invokeRestart("muffleWarning")
      }
    )
  },
  error = function(e){
    message(e$message)
    return(NA)
  })
}

cor.test.wErrorHandling <- function(x,y,...){
  tryCatch({
    withCallingHandlers(
      {
        cor.test(x,y,...)
      },
      warning = function(w){
        invokeRestart("muffleWarning")
      }
    )
  },
  error = function(e){
    message(e$message)
    estimateNumber=NA
    names(estimateNumber)=NA
    return(list(estimate=estimateNumber,p.value=NA))
  })
}

wilcox.test.p <- function(x,y,...){
  result <- wilcox.test.wErrorHandling(x,y,...)
  if(any(is.na(result))){return(NA)}else{return(result$p.value)}}

getDifferentialAlterationByStageData = function(diffExpData, CNAdata, sampleData, thisCohort, thisSubtype, p.adjust.method = "fdr", stage1="IBC", stage2="LVI"){
  
  #Complete the Diff Exp data with gene information
  #################################################
  #First, we get all by geneID and remove duplicates
  diffExpEntrez <- inner_join(diffExpData,
                              grch38 %>% filter(chr %in% validCHRs), #we filter out genes in alternative scaffolds since they are not useful for this
                              by = join_by(entrez==entrez),
                              na_matches = "never", 
                              keep = T, #Keep both keys with the suffix below
                              suffix = c(".diffExp",".DB"),
                              multiple = "all") %>%
    group_by(entrez.diffExp) %>% #If some had multiple matches, we sort the groups per gene size and keep the top (largest)
    arrange(desc(end-start),.by_group = T) %>%
    slice(1)
  
  #Then, we get them by name, remove duplicates, and can complete the list using names, only for those names not found yet
  diffExpSymbol <- inner_join(diffExpData,
                              grch38 %>% filter(chr %in% validCHRs), #we filter out genes in alternative scaffolds since they are not useful for this
                              by = join_by(symbol==symbol),
                              na_matches = "never", 
                              keep = T, #Keep both keys with the suffix below
                              suffix = c(".diffExp",".DB"),
                              multiple = "all") %>%
    group_by(symbol.diffExp) %>% #If some had multiple matches, we sort the groups per gene size and keep the top (largest)
    arrange(desc(end-start),.by_group = T) %>%
    slice(1)
  
  #Concatenate the ones found by name that were not found yet
  diffExpFinal <- rbind(diffExpEntrez,
                        anti_join(diffExpSymbol,diffExpEntrez,na_matches = "never",by=join_by(symbol.diffExp)))
  
  allResults <- bind_rows(mclapply(c(1:nrow(diffExpFinal)),FUN = function(iGene){
    thisGene <- diffExpFinal[iGene,]
    CNAdata %>% 
      filter(chromosome==thisGene$chr & end>thisGene$start & start<thisGene$end) %>%
      mutate(overlap=calcOverlap(thisGene$start,thisGene$end,start,end),symbol.diffExp=thisGene$symbol.diffExp) %>%
      group_by(sample) %>%
      summarize(across(-c("overlap","logR","chromosome","start","end"),first),WeightedMeanLogR = weighted.mean(logR,w = overlap/sum(overlap)))
  },mc.cores=mc_cores))
  
  allResultsMerged <- left_join(allResults,diffExpFinal,by=join_by(symbol.diffExp))
  allResultsMerged <- left_join(allResultsMerged,sampleData%>%mutate(patient=as.character(patient)),by=join_by(sample==fullLabel,patient==patient),keep = F)
   
  #Paired samples, selecting using the variance
  #First we select the samples we will use
  #We are using the data merged with the gene CNA information to select samples taking into account their missing data
  selectedSamples <- allResultsMerged %>%
    #Considering missing data
    select(c(sample,patient,Stage,Cohort,LVI_subtype,observedVariance,symbol.diffExp)) %>%
    count(across(c(sample,patient,Stage,Cohort,LVI_subtype,observedVariance)),sort = T,name = "nGenesWithData") %>%
    #Filtering patients and samples of interest
    filter(Cohort==thisCohort) %>%
    {if (!is.null(thisSubtype)) filter(.,LVI_subtype==thisSubtype) else . } %>%
    filter(Stage%in%c(stage1,stage2)) %>% #Only considering 2 stages
    #Grouping by patient
    group_by(patient) %>% 
    filter(n_distinct(Stage) == 2) %>% #Only patients with samples from 2 stages #Without this we can not do them paired but we don't have dependent samples within group
    group_by(patient,Stage) %>%
    #slice(sample(n())) %>% #Reorder tibble so that we are selecting the samples at random if not using variance
    #arrange(desc(nGenesWithData),.by_group = T) %>%
    #But then sort by data completeness 
    arrange(desc(nGenesWithData),observedVariance,.by_group = T) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  pairedResults <- allResultsMerged %>%
    semi_join(selectedSamples, by = "sample") %>% 
    group_by(symbol.diffExp) %>%
    summarize(
      meanStage1 = mean(WeightedMeanLogR[Stage == stage1]),
      meanStage2 = mean(WeightedMeanLogR[Stage == stage2]),
      across(c(entrez.diffExp,category,log2FoldChange),unique),
      wilcox.pval = wilcox.test.p(WeightedMeanLogR[Stage == stage1],WeightedMeanLogR[Stage == stage2],paired = T)) %>%
    mutate(log2FoldChangeCNA = meanStage2-meanStage1, wilcox.adjusted.pval = p.adjust(wilcox.pval,method=p.adjust.method))
  
  #Not paired only using one samples per group and patient
  selectedSamples <- allResultsMerged %>%
    #Considering missing data
    select(c(sample,patient,Stage,Cohort,LVI_subtype,observedVariance,symbol.diffExp)) %>%
    count(across(c(sample,patient,Stage,Cohort,LVI_subtype,observedVariance)),sort = T,name = "nGenesWithData") %>%
    #Filtering patients and samples of interest
    filter(Cohort==thisCohort) %>%
    {if (!is.null(thisSubtype)) filter(.,LVI_subtype==thisSubtype) else . } %>%
    filter(Stage%in%c(stage1,stage2)) %>% #Only considering 2 stages
    #Grouping by patient
    group_by(patient) %>% 
    group_by(patient,Stage) %>%
    #slice(sample(n())) %>% #Reorder tibble so that we are selecting the samples at random if not using variance
    #arrange(desc(nGenesWithData),.by_group = T) %>%
    #But then sort by data completeness 
    arrange(desc(nGenesWithData),observedVariance,.by_group = T) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  notPairedResults <- allResultsMerged %>%
    semi_join(selectedSamples, by = "sample") %>% 
    group_by(symbol.diffExp) %>%
    summarize(
      meanStage1 = mean(WeightedMeanLogR[Stage == stage1]),
      meanStage2 = mean(WeightedMeanLogR[Stage == stage2]),
      across(c(entrez.diffExp,category,log2FoldChange),unique),
      wilcox.pval = wilcox.test.p(WeightedMeanLogR[Stage == stage1],WeightedMeanLogR[Stage == stage2],paired = F)) %>%
    mutate(log2FoldChangeCNA = meanStage2-meanStage1, wilcox.adjusted.pval = p.adjust(wilcox.pval,method=p.adjust.method))
  
  return(list(notPaired=notPairedResults,paired=pairedResults))
}

getDifferentialAlterationBySubtypeData = function(diffExpData, CNAdata, sampleData, thisCohort, thisStage, p.adjust.method = "fdr", subtype1="ID", subtype2="Prol"){
  
  #Complete the Diff Exp data with gene information
  #################################################
  #First, we get all by geneID and remove duplicates
  diffExpEntrez <- inner_join(diffExpData,
                              grch38 %>% filter(chr %in% validCHRs), #we filter out genes in alternative scaffolds since they are not useful for this
                              by = join_by(entrez==entrez),
                              na_matches = "never", 
                              keep = T, #Keep both keys with the suffix below
                              suffix = c(".diffExp",".DB"),
                              multiple = "all") %>%
    group_by(entrez.diffExp) %>% #If some had multiple matches, we sort the groups per gene size and keep the top (largest)
    arrange(desc(end-start),.by_group = T) %>%
    slice(1)
  
  #Then, we get them by name, remove duplicates, and can complete the list using names, only for those names not found yet
  diffExpSymbol <- inner_join(diffExpData,
                              grch38 %>% filter(chr %in% validCHRs), #we filter out genes in alternative scaffolds since they are not useful for this
                              by = join_by(symbol==symbol),
                              na_matches = "never", 
                              keep = T, #Keep both keys with the suffix below
                              suffix = c(".diffExp",".DB"),
                              multiple = "all") %>%
    group_by(symbol.diffExp) %>% #If some had multiple matches, we sort the groups per gene size and keep the top (largest)
    arrange(desc(end-start),.by_group = T) %>%
    slice(1)
  
  #Concatenate the ones found by name that were not found yet
  diffExpFinal <- rbind(diffExpEntrez,
                        anti_join(diffExpSymbol,diffExpEntrez,na_matches = "never",by=join_by(symbol.diffExp)))
  
  allResults <- bind_rows(mclapply(c(1:nrow(diffExpFinal)),FUN = function(iGene){
    thisGene <- diffExpFinal[iGene,]
    CNAdata %>% 
      filter(chromosome==thisGene$chr & end>thisGene$start & start<thisGene$end) %>%
      mutate(overlap=calcOverlap(thisGene$start,thisGene$end,start,end),symbol.diffExp=thisGene$symbol.diffExp) %>%
      group_by(sample) %>%
      summarize(across(-c("overlap","logR","chromosome","start","end"),first),WeightedMeanLogR = weighted.mean(logR,w = overlap/sum(overlap)))
  },mc.cores=mc_cores))
  
  allResultsMerged <- left_join(allResults,diffExpFinal,by=join_by(symbol.diffExp))
  allResultsMerged <- left_join(allResultsMerged,sampleData%>%mutate(patient=as.character(patient)),by=join_by(sample==fullLabel,patient==patient),keep = F)
   
  #Not paired only using one sample per group and patient
  selectedSamples <- allResultsMerged %>%
    #Considering missing data
    select(c(sample,patient,Stage,Cohort,LVI_subtype,observedVariance,symbol.diffExp)) %>%
    count(across(c(sample,patient,Stage,Cohort,LVI_subtype,observedVariance)),sort = T,name = "nGenesWithData") %>%
    #Filtering patients and samples of interest
    filter(Cohort==thisCohort) %>%
    {if (!is.null(thisStage)) filter(.,Stage==thisStage) else .} %>%
    filter(LVI_subtype%in%c(subtype1,subtype2)) %>% #Only considering 2 subtypes
    #Grouping by patient
    group_by(patient) %>% 
    group_by(patient,LVI_subtype) %>%
    #slice(sample(n())) %>% #Reorder tibble so that we are selecting the samples at random if not using variance
    #arrange(desc(nGenesWithData),.by_group = T) %>%
    #But then sort by data completeness 
    arrange(desc(nGenesWithData),observedVariance,.by_group = T) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  notPairedResults <- allResultsMerged %>%
    semi_join(selectedSamples, by = "sample") %>% 
    group_by(symbol.diffExp) %>%
    summarize(
      meanSubtype1 = mean(WeightedMeanLogR[LVI_subtype == subtype1]),
      meanSubtype2 = mean(WeightedMeanLogR[LVI_subtype == subtype2]),
      across(c(entrez.diffExp,category,log2FoldChange),unique),
      wilcox.pval = wilcox.test.p(WeightedMeanLogR[LVI_subtype == subtype1],WeightedMeanLogR[LVI_subtype == subtype2],paired = F)) %>%
    mutate(log2FoldChangeCNA = meanSubtype2-meanSubtype1, wilcox.adjusted.pval = p.adjust(wilcox.pval,method=p.adjust.method))
  
  return(list(notPaired=notPairedResults))
}


printResults <- function(theseResults, alpha=0.05, maxRows=1000, pvalcolname="wilcox.adjusted.pval", pvalcolnameCorr="wilcox.pval", printPairings = c("paired","notPaired"), cor.method="kendall"){
  for (condition in names(theseResults)) {
    for (pairing in printPairings){
      if(nrow(theseResults[[condition]][[pairing]]%>%filter(!!sym(pvalcolname)<alpha))>0){
        print(sprintf("The following differentially-expressed genes for %s using the %s strategy were also differentially copy-number altered",condition,pairing))
        theseResults[[condition]][[pairing]]%>%filter(!!sym(pvalcolname)<alpha) %>% arrange(!!sym(pvalcolname)) %>% print(n=maxRows)
      } else {
        print(sprintf("None of the differentially-expressed genes for %s using the %s strategy were differentially copy-number altered",condition,pairing))
      }
      allCorrelation <- cor.test.wErrorHandling(theseResults[[condition]][[pairing]] %>% pull(log2FoldChange),
                                 theseResults[[condition]][[pairing]] %>% pull(log2FoldChangeCNA),method = cor.method)
      sigCorrelation <- cor.test.wErrorHandling(theseResults[[condition]][[pairing]] %>% filter(!!sym(pvalcolnameCorr)<alpha) %>% pull(log2FoldChange),
                                 theseResults[[condition]][[pairing]] %>% filter(!!sym(pvalcolnameCorr)<alpha) %>% pull(log2FoldChangeCNA),method = cor.method)
      writeLines(sprintf("logR correlation: %s = %.2f, p = %e; logR correlation for %s < %.2f: %s = %.2f, p = %e \n",
                         names(allCorrelation$estimate),
                         allCorrelation$estimate,
                         allCorrelation$p.value, 
                         pvalcolnameCorr, 
                         alpha, 
                         names(sigCorrelation$estimate),
                         sigCorrelation$estimate,
                         sigCorrelation$p.value))
    }
  }
}

#Configuration
configFile <- paste(sep="/",Sys.getenv("lviProjectDNAScripts"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the lviProjectDNAScripts environment variable before running this script")
}
source(configFile)
inDir <- acnaDir
varFile <- paste(sep="/",acnaDir,"variance.tsv")

#Files
pattern <- "^(.*)[:.:]([^.]*)[:.:]best[:.:]ACNAs[:.:]c(.*)[:.:]p(.*)[:.:]tsv$" ##[::] is R's way of escaping something, similar to \ or \\ in most regex
#Samples discarded due to problems
#These normals clearly do not represent a normal sample, with plenty alterations very similar to the Met in the same patient
badSamples <- tibble(patient=c("344","740"),sample=c("NORMAL_2","NORMAL_1"))

#CNA data
##################

#IO of previous results with absolute and relative copy number calls using QDNAseq + Rascal
fileList <- dir(inDir, pattern, full.names = TRUE)
allCalls <- bind_rows(lapply(fileList,function(x,pattern){
  theData <- read_delim(x,show_col_types = F);
  theData %>% mutate(patient=gsub(basename(x),pattern=pattern,replacement="\\1"),
                     sample=gsub(basename(x),pattern=pattern,replacement="\\2"),
                     cellularity=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\3")),
                     ploidy=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\4")))},pattern=pattern))

CNAdata <- anti_join(allCalls,badSamples,by = c("patient","sample")) %>% #Filter bad samples
  mutate(logR = log2(ifelse(absolute_copy_number<0,0.1,absolute_copy_number)/ploidy)) %>% #logR calculated using absolute copy number and ploidy
  rename(patientSample = sample) %>%
  mutate(sample = paste(sep=".",patient,patientSample))
  #mutate(logR = ifelse(round(ploidy) == round(absolute_copy_number),0,logR)) #logR of regions estimated to be not altered fixed to 0 (not a good idea here?)

validCHRs <- CNAdata %>% pull(chromosome) %>% unique

#Patient data
##################
sampleData <- read_tsv(sampleDataFile, show_col_types = F)
varData <- read_delim(varFile, show_col_types = F)
sampleData <- full_join(sampleData,varData,by = c("patient","label"))

#Genes of interest
##################

#Stages
#######
transitionFrom <- c("NORMAL","DCIS","IBC","LVI")
transitionTo <- c("DCIS","IBC","LVI","MET")
transitions <- unite(tibble(transitionFrom,transitionTo),"transition",sep="to") %>% pull(transition)

discoveryData <- lapply(1:length(transitions),function(x){
  read_xlsx(paste(sep="/",diffExpDir,diffExpStageDisc),sheet = x) %>% 
    mutate(entrez = as.numeric(entrez), category = ifelse(log2FoldChange>0,"Up","Down")) %>%
    rename(symbol=Symbol)})

names(discoveryData) <- transitions

allResultsDiscoveryStages <- lapply(1:length(discoveryData),function(iDataset){
  return(getDifferentialAlterationByStageData(discoveryData[[iDataset]],CNAdata,sampleData,thisCohort = "Disc",thisSubtype = NULL,stage1 = transitionFrom[iDataset],stage2 = transitionTo[iDataset]))})

names(allResultsDiscoveryStages) <- transitions

printResults(allResultsDiscoveryStages,printPairings = "paired")
#printResults(allResultsDiscoveryStages,pvalcolname = "wilcox.pval")

validationData <- lapply(1:length(transitions),function(x){
  read_xlsx(paste(sep="/",diffExpDir,diffExpStageVal), sheet = x) %>% 
    mutate(entrez = as.numeric(entrez), category = ifelse(log2FoldChange>0,"Up","Down")) %>%
    rename(symbol=Symbol)})

names(validationData) <- transitions

allResultsValidationStages <- lapply(1:length(validationData),function(iDataset){
  return(getDifferentialAlterationByStageData(validationData[[iDataset]],CNAdata,sampleData,thisCohort = "Val",thisSubtype = NULL,stage1 = transitionFrom[iDataset],stage2 = transitionTo[iDataset]))})

names(allResultsValidationStages) <- transitions

printResults(allResultsValidationStages,printPairings = "paired")
#printResults(allResultsValidationStages,pvalcolname = "wilcox.pval")

#LVI across subtype
###################

LVIAcrossDisc <- read_delim(paste(sep="/",diffExpDir,diffExpAcrossDisc), show_col_types = F) %>%
  mutate(entrez = as.numeric(entrez), category = ifelse(log2FoldChange>0,"Up","Down")) %>%
  rename(symbol = "...1")

resultsLVIAcrossDisc <- getDifferentialAlterationBySubtypeData(LVIAcrossDisc,CNAdata,sampleData,thisCohort = "Disc",thisStage = "LVI",subtype1 = "ID",subtype2 = "Prol")

LVIAcrossVal <- read_delim(paste(sep="/",diffExpDir,diffExpAcrossVal), show_col_types = F) %>%
  mutate(entrez = as.numeric(entrez), category = ifelse(log2FoldChange>0,"Up","Down")) %>%
  rename(symbol = "...1")

resultsLVIAcrossVal <- getDifferentialAlterationBySubtypeData(LVIAcrossVal,CNAdata,sampleData,thisCohort = "Val",thisStage = "LVI",subtype1 = "Prol",subtype2 = "ID")

resultsLVIAcross <- list(Disc=resultsLVIAcrossDisc,Val=resultsLVIAcrossVal)

printResults(resultsLVIAcross,printPairings = "notPaired")


#LVI vs IBC by subtype
######################
DiscIDDiffExpDown <- read_xlsx(paste(sep="/",diffExpDir,diffExpIBCLVIDiscID), sheet = 1)
DiscIDDiffExpUp <- read_xlsx(paste(sep="/",diffExpDir,diffExpIBCLVIDiscID), sheet = 2)
DiscIDDiffExp <- rbind(DiscIDDiffExpDown %>% mutate(category="Down"),
      DiscIDDiffExpUp %>% mutate(category="Up")) %>% mutate(entrez = as.numeric(entrez)) %>%
  rename(symbol=Symbol)

DiscProlDiffExp <- read_xlsx(paste(sep="/",diffExpDir,diffExpIBCLVIDiscNID), sheet = 1)
DiscProlDiffExp <- DiscProlDiffExp %>% 
  mutate(category=ifelse(log2FoldChange>0,"Up","Down")) %>%
  mutate(entrez = as.numeric(entrez)) %>%
  rename(symbol = Symbol)

ValIDDiffExp <- read_csv(paste(sep="/",diffExpDir,diffExpIBCLVIValID), show_col_types = F) %>%
  rename(symbol = "...1") %>% 
  #Direction inverted to make it the same as in the Discovery cohort
  mutate(category=ifelse(log2FoldChange<0,"Up","Down"), log2FoldChange = log2FoldChange*-1)

ValProlDiffExp <- read_csv(paste(sep="/",diffExpDir,diffExpIBCLVIValNID), show_col_types = F) %>%
  rename(symbol = "...1") %>% 
  #Direction inverted to make it the same as in the Discovery cohort
  mutate(category=ifelse(log2FoldChange<0,"Up","Down"), log2FoldChange = log2FoldChange*-1)

datasetCohorts <- c("Disc","Disc","Val","Val")
datasetSubtypes <- c("ID","Prol","ID","Prol")
datasets <- list(DiscID = DiscIDDiffExp,DiscProl = DiscProlDiffExp,ValID = ValIDDiffExp,ValProl = ValProlDiffExp)


allResultsSubtype <- lapply(1:length(datasets),function(iDataset){
  return(getDifferentialAlterationByStageData(datasets[[iDataset]],CNAdata,sampleData,datasetCohorts[iDataset],datasetSubtypes[iDataset]))})
names(allResultsSubtype) <- names(datasets)

printResults(allResultsSubtype,printPairings = "paired")
#printResults(allResultsSubtype,pvalcolname = "wilcox.pval")

